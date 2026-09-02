module Force
  use allocation
  use Geometry

  contains


  subroutine Force_Calculation

    implicit none

    integer :: prev_idx, next_idx
    real*8 :: grad_perimeter_X, grad_perimeter_Y, grad_area_X, grad_area_Y
    real*8 :: rann
    ! OPTIMIZATION (log.txt): edge_dx/edge_dy/edge_len(k) cache the vector and
    ! length of the edge from vertex k to vertex next(k), computed once per
    ! cell instead of twice (each edge was previously recomputed once as
    ! vertex jc's "next" edge and again as vertex jc+1's "prev" edge, halving
    ! the sqrt calls in this hot loop -- run every cell, every timestep).
    real*8 :: edge_dx(inn_dim1), edge_dy(inn_dim1), edge_len(inn_dim1)


    fxx = 0.0d0
    fyy = 0.0d0

    ! BUGFIX (log.txt): beta/gamm nondimensionalization moved to read_input
    ! (allocation.f90), done once instead of every call -- see log.txt for why
    ! repeating it here compounded beta/gamm every timestep.

    do ic = 1, Nc !Lx*Ly
      
      nn = num(ic)

      vx = v(1,inn(1:nn,ic))
      vy = v(2,inn(1:nn,ic))
      

      if(if_RhoROCK)then
        if(if_coupling_noise)then
          call random_number(rann)
          rann = rann - 0.5d0
          ! BUGFIX (log.txt): this branch omitted the /(lambda*Ao) factor that
          ! the else-branch below applies, silently changing beta's scale by
          ! a factor of (lambda*Ao) whenever coupling noise is enabled.
          beta = (Myosin_Coupling_Strength * Myosin(ic) + coupling_noise_strength * rann)/(lambda*Ao)
        else
          beta = Myosin_Coupling_Strength * Myosin(ic)/(lambda*Ao)
        end if
      end if

      if(if_active_contractility)then
        beta = Myosin(ic) / (lambda * Ao)
      end if



      call CalculateArea(vx,vy,nn,area)
      call CalculatePerimeter(vx,vy,nn,perimeter)
      ! NOTE (log.txt): CalculateArea returns the signed shoelace area; only
      ! the magnitude is meaningful here (2D area has no physical sign), so
      ! area is intentionally take as abs(area). Confirmed with user: mesh
      ! winding is a fixed convention in this codebase, no per-cell sign
      ! handling needed in the gradient below.
      area = abs(area)

      ! edge_dx/dy/len(k) = vector/length of the edge from vertex k to
      ! vertex k's next neighbour (wraparound at nn).
      do jc = 1, nn
        next_idx = jc + 1
        if(jc == nn) next_idx = 1
        edge_dx(jc) = vx(next_idx) - vx(jc)
        edge_dy(jc) = vy(next_idx) - vy(jc)
        edge_len(jc) = dsqrt(edge_dx(jc)**2 + edge_dy(jc)**2)
      end do

      do jc = 1,num(ic)

       next_idx = jc + 1
       prev_idx = jc - 1
       if(jc == num(ic))then
         next_idx = 1
       elseif(jc == 1)then
         prev_idx = num(ic)
       end if

       if(prev_idx.eq.0)then
         write(*,*)'OOPS!!!', ic, jc, num(ic), inn(1:num(ic), ic)
       end if

        grad_area_X = 0.5d0 * (vy(prev_idx) - vy(next_idx))
        grad_area_Y = 0.5d0 * (vx(next_idx) - vx(prev_idx))

        ! dx/dy/len_d for the (prev_idx -> jc) edge = edge_dx/dy/len(prev_idx);
        ! for the (jc -> next_idx) edge = -edge_dx/dy(jc), edge_len(jc).
        grad_perimeter_X = edge_dx(prev_idx)/edge_len(prev_idx) - edge_dx(jc)/edge_len(jc)
        grad_perimeter_Y = edge_dy(prev_idx)/edge_len(prev_idx) - edge_dy(jc)/edge_len(jc)


        fxx(inn(jc,ic)) = fxx(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_X &
          - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_X &
          - gamm * grad_perimeter_X

        fyy(inn(jc,ic)) = fyy(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_Y  &
          - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_Y &
          - gamm * grad_perimeter_Y

      end do

   end do

   fxx = fxx 
   fyy = fyy 

  end subroutine Force_Calculation


  subroutine Motile_Force_Calculation

    implicit none

    integer :: nv
    ! OPTIMIZATION (log.txt): fixed-size automatic (stack) arrays sized to
    ! the reserved capacity v_dim2 -- only the first nv elements are ever
    ! touched. Declared once per call, no heap allocate/deallocate.
    real*8 :: rann_x(v_dim2), rann_y(v_dim2), sigma_mot(v_dim2)

    ! OPTIMIZATION (log.txt): mot/fxx_ran/fyy_ran are sized v_dim2, the
    ! RESERVED vertex capacity for future cell-division growth -- not the
    ! number of vertices actually in use right now. Vertices are only ever
    ! appended (Proliferation_Core writes new ones at maxval(inn)+1/+2,
    ! never compacted or reused), so maxval(inn) is a correct, cheap
    ! (O(inn_dim1*inn_dim2) integer scan) upper bound on which vertex
    ! indices are actually referenced by any live cell. Slots beyond it are
    ! never read by anything -- CalculateArea/CalculatePerimeter/
    ! Force_Calculation only ever iterate via inn(1:nn,ic).
    ! NOTE: with if_cell_division off (or growth far short of the reserved
    ! capacity), this trims only a small, fixed number of always-unused
    ! slots -- the much bigger win below is eliminating ~2*nv individual
    ! random_number() calls (measured as the dominant cost in profiling of
    ! a 10^7-step run) in favor of 2 whole-array calls, which removes the
    ! per-call subroutine overhead of tens of thousands of scalar calls per
    ! timestep, plus lets the compiler vectorize the arithmetic below.
    ! NOTE ON REPRODUCIBILITY: unlike the sqrt-caching/early-exit
    ! optimizations elsewhere in this codebase, this one is NOT guaranteed
    ! to reproduce the exact same trajectory as the old scalar-loop version
    ! for a given (unseeded) run -- both trimming to nv<v_dim2 and batching
    ! the random_number calls change how many/which calls consume the
    ! global RNG stream, which shifts everything drawn afterward (this
    ! codebase has no random_seed anywhere, so exact run-to-run
    ! reproducibility was never guaranteed regardless -- see log.txt). The
    ! physics is unaffected: every live vertex still gets an independent,
    ! identically-distributed noise term.
    nv = maxval(inn)

    if(if_motility_decay)then
     ! mot = mot0*exp(-it*dt/motility_decay_timeScale)

     mot(1:nv) = mot(1:nv) - dt * mot(1:nv)/motility_decay_timeScale
    end if

    call random_number(rann_x(1:nv))
    call random_number(rann_y(1:nv))

    ! sqrt(2*mot*eta) is identical for both components -- was computed
    ! twice per vertex (once for fxx_ran, once for fyy_ran); now once.
    sigma_mot(1:nv) = sqrt(2.0d0 * mot(1:nv) * eta)

    fxx_ran(1:nv) = sigma_mot(1:nv) * (2.0d0 * rann_x(1:nv) - 1.0d0)
    fyy_ran(1:nv) = sigma_mot(1:nv) * (2.0d0 * rann_y(1:nv) - 1.0d0)

  end subroutine Motile_Force_Calculation

  subroutine ABP_Force_Calculation
    implicit none

    real*8 :: rann
    
    if (if_ABP)then

      do ic = 1, Nc

        ! NOTE (log.txt): one rann draw per cell here is intentional -- ABP
        ! treats each cell as a single self-propelled particle with one
        ! rotational-noise process, applied uniformly to that particle's own
        ! vertices. This is a distinct concept from Polar_Motile_Force_
        ! Calculation (a per-vertex force that must accumulate contributions
        ! from every cell touching a shared vertex) -- the two should not be
        ! reconciled to the same pattern. Reverted the earlier per-vertex
        ! noise change on that basis.
        call gaussian_random(rann, 0.0d0, 1.0d0)  ! number, mu, sigma

        do jc = 1, num(ic)

          fxx_ABP(inn(jc,ic)) = vo*dcos(theta_ABP(inn(jc,ic)))
          fyy_ABP(inn(jc,ic)) = vo*dsin(theta_ABP(inn(jc,ic)))
          rot_noise(inn(jc,ic)) = sqrt(2.0d0 * Dr) * rann

        end do

      end do

    else
      fxx_ABP = 0.0d0
      fyy_ABP = 0.0d0
      theta_ABP = 0.0d0
      rot_noise = 0.0d0
    end if


  end subroutine ABP_Force_Calculation

  subroutine Polar_Motile_Force_Calculation
     implicit none
     real*8 :: rann_x, rann_y

     if(if_polar_motility)then

       ! BUGFIX (log.txt): fxx_Polar/fyy_Polar were assigned with '=' instead
       ! of accumulated with '+', unlike every other per-vertex force in this
       ! module -- a vertex shared between cells simply got overwritten by
       ! whichever cell was processed last, discarding the others' contribution.
       fxx_Polar = 0.0d0
       fyy_Polar = 0.0d0

       do ic = 1, Nc
         call random_number(rann_x)
         call random_number(rann_y)

         do jc = 1, num(ic)
           fxx_Polar(inn(jc,ic)) = fxx_Polar(inn(jc,ic)) + & 
             polar_motility_strength * (rann_x - 0.50d0)
           fyy_Polar(inn(jc,ic)) = fyy_Polar(inn(jc,ic)) + & 
             polar_motility_strength * (rann_y - 0.50d0)
         end do
       end do


     else
       fxx_Polar = 0.0d0
       fyy_Polar = 0.0d0
     end if


  end subroutine Polar_Motile_Force_Calculation

  subroutine Give_Motility_Gradient

    implicit none


    integer :: ip, idx
    real(8) :: vy, lowpatch, highpatch, meanpatch
    real(8) :: patchwidth
    integer :: verryCount(size(v, 2))
    ! BUGFIX (log.txt): verry stores a real*8 y-coordinate (v(2,...)) and must
    ! not be declared integer -- it was silently truncating the fractional
    ! part of every vertex y-position before it fed into the patch binning.
    real(8) :: verry(size(v, 2))
    real(8), allocatable :: veryyN(:)
    integer :: numver
    integer, allocatable :: valid_indices(:)





    ! Initialize arrays
    verryCount = 0
    verry = 0

    ! Loop over the cells
    do ic = 1, Nc !Lx * Ly
        do jc = 1, num(ic)
            verry(inn(jc, ic)) = v(2, inn(jc, ic))
            verryCount(inn(jc, ic)) = 1
        end do
    end do


    ! Find valid indices where verryCount is non-zero
    !valid_indices = 0
    idx = 1
    do ip = 1, size(verryCount)
        if (verryCount(ip) .ne. 0) then
            !valid_indices(idx) = ip
            idx = idx + 1
        end if
    end do

    allocate(valid_indices(idx))
    valid_indices = 0
    
    idx = 1
    do ip = 1, size(verryCount)
        if (verryCount(ip) .ne. 0) then
            valid_indices(idx) = ip
            idx = idx + 1
        end if
    end do


    ! Extract valid values for veryyN based on valid_indices
    numver = idx - 1
    if (numver > 0) then
        allocate(veryyN(numver))
        do ip = 1, numver
            veryyN(ip) = verry(valid_indices(ip))
        end do
    end if

    patchwidth = 1.0d0
    lowpatch = minval(veryyN) - patchwidth / 2.0d0
    highpatch = lowpatch + patchwidth / 2.0d0
    meanpatch = (highpatch + lowpatch) / 2.0d0

    ! Initialize mot array
    mot = 0.0

    ! Loop over patches
    do ip = 1, numver
        do ic = 1, Nc !Lx * Ly
            do jc = 1, num(ic)
                vy = v(2, inn(jc, ic))
                if (lowpatch <= vy .and. vy < highpatch) then
                    mot(inn(jc, ic)) = etas_max * exp(-vy / (mot_Lc*Ly))
                end if
            end do
        end do

        lowpatch = highpatch
        highpatch = lowpatch + patchwidth
        meanpatch = (highpatch + lowpatch) / 2.0d0
    end do

    deallocate(veryyN)
    deallocate(valid_indices)

end subroutine Give_Motility_Gradient


subroutine Give_Motility_Hotspot

  implicit none
  integer :: ic, jc, ip, jp

!  integer, parameter :: number_of_hotspot = 4
  integer, dimension(:), allocatable :: hotspot_location  ! Cell index
  real*8, dimension(number_of_hotspot) :: xCM, yCM

!  real*8 :: sigma_hotspot
  real*8 :: xij, yij, rij

  print*,'number of hotspot, sigma_hotspot', number_of_hotspot, sigma_hotspot

  if(number_of_hotspot.eq.4)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ 12, Ly-12, Lx*Ly-Ly+1+12,  Lx*Ly-12 /)

  elseif(number_of_hotspot.eq.2)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ 12, Lx*Ly-12/)

  elseif(number_of_hotspot.eq.1)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ int(dble(Ly)/2.0d0)/)

  else
    ! BUGFIX (log.txt): without this branch, an unsupported number_of_hotspot
    ! left hotspot_location unallocated, and it is dereferenced unconditionally
    ! below (ip = hotspot_location(jp)) -- a crash. Fail loudly instead.
    write(*,*) 'ERROR: number_of_hotspot must be 1, 2, or 4. Got:', number_of_hotspot
    stop 1

  end if

!  print*, hotspot_location
!  print*, xCM
!  print*, yCM

  do jp = 1, number_of_hotspot
    ip = hotspot_location(jp)
    xCM(jp) = sum(v(1,inn(1:num(ip),ip)))/dble(num(ip))
    yCM(jp) = sum(v(2,inn(1:num(ip),ip)))/dble(num(ip))
  end do
  
 ! sigma_hotspot = 5.0d0


  do ic  = 1, Nc !Lx*Ly
    do jc = 1, num(ic) 
      
      do ip = 1, number_of_hotspot

        
        xij = v(1,inn(jc,ic)) - xCM(ip)
        yij = v(2,inn(jc,ic)) - yCM(ip)
        rij = xij**2 + yij**2
    
        mot(inn(jc,ic)) = mot(inn(jc,ic)) + & 
          etas_max * exp(-(rij)/sigma_hotspot**2) / 3.0d0
      end do

    end do
  end do

end subroutine Give_Motility_Hotspot

subroutine Apply_Fixed_Boundary

  implicit none

  integer :: im

  if(sum(Total_T2_count(1:it)).gt.0)then
     call Find_boundary_dynamic
   else 
     call Get_Boundary_info
   end if

  do ii = 1,size(boundary)
     im = boundary(ii)
     nn = num(im)
     fxx(inn(1:nn, im)) = 0.0d0
     fxx_ran(inn(1:nn, im)) = 0.0d0
     fyy(inn(1:nn, im)) = 0.0d0
     fyy_ran(inn(1:nn, im)) = 0.0d0
     fxx_ABP(inn(1:nn, im)) = 0.0d0
     fyy_ABP(inn(1:nn, im)) = 0.0d0
     fxx_Polar(inn(1:nn, im)) = 0.0d0
     fyy_Polar(inn(1:nn, im)) = 0.0d0

   end do


end subroutine Apply_Fixed_Boundary

subroutine Apply_bottom_border_Fixed
  implicit none

!  call Find_boundary_dynamic

!  print*, 'Bottom border fixed'

  fxx(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_ran(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_ran(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_ABP(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_ABP(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_Polar(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_Polar(bottom_border(1:bottom_border_count)) = 0.0d0


end subroutine  Apply_bottom_border_Fixed

subroutine Apply_top_border_Fixed
  implicit none

!  print*, 'Top border fixed'

  fxx(top_border(1:top_border_count)) = 0.0d0
  fxx_ran(top_border(1:top_border_count)) = 0.0d0
  fyy(top_border(1:top_border_count)) = 0.0d0
  fyy_ran(top_border(1:top_border_count)) = 0.0d0
  fxx_ABP(top_border(1:top_border_count)) = 0.0d0
  fyy_ABP(top_border(1:top_border_count)) = 0.0d0
  fxx_Polar(top_border(1:top_border_count)) = 0.0d0
  fyy_Polar(top_border(1:top_border_count)) = 0.0d0


end subroutine Apply_top_border_Fixed

subroutine Apply_perturbation
  implicit none
  real*8 :: comb


  if(if_sin_perturb)then
    if(it.ge.sin_perturb_when)then
       if(if_dirac_comb)then

          if(mod((it-sin_perturb_when), & 
            comb_onPeriod+comb_offPeriod) &
            < comb_onPeriod)then

            comb = 1.0d0
          else 
            comb = 0.0d0
          end if

          v(1,:) = v(1,:) + dt * sin_perturb_strength * &
            sin(2*pi*sin_perturb_waveNumber* v(2,:)/Ly) * comb


       end if
    end if
  end if




end subroutine Apply_perturbation


subroutine Apply_Limb_Force
! Limb force at four spots near the boundary, towards outward.  


  implicit none 

  integer :: ic, ip

  integer, dimension(4) :: applyforce_at


    applyforce_at = (/ 4, Ly-4, Lx*Ly-Ly+1+4,  Lx*Ly-4 /)

    do ip =  1, size(applyforce_at)
      ic = applyforce_at(ip)
      nn = num(ic)

      ! BUGFIX (log.txt): every other force in the code converts to a
      ! displacement via dt*force/eta (see vertexmain.f90); this was missing
      ! the /eta mobility factor, making its effect inconsistent by a factor
      ! of eta relative to every other named force.
      if(ic .lt. int(Lx*Ly/2))then   ! Cheap way
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (-1.0d0) / eta
      else
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (1.0d0) / eta
      end if

    end do


end subroutine Apply_Limb_Force

subroutine gaussian_random(rnumber, mu, sigma)
  implicit none
  real(8), intent(out) :: rnumber
  real(8), intent(in), optional :: mu, sigma

  real(8) :: u1, u2, r, theta, z
  real(8), parameter :: pi = 4.0d0 * datan(1.0d0)
  real(8) :: mean, std

  ! Set defaults
  mean = 0.0d0
  std  = 1.0d0
  if (present(mu)) mean = mu
  if (present(sigma)) std = sigma

  ! Generate two uniform random numbers in (0,1)
  call random_number(u1)
  call random_number(u2)

  ! Ensure u1 is not zero (log(0) issue)
  if (u1 < 1.0d-12) u1 = 1.0d-12

  ! Box–Muller transform
  r = dsqrt(-2.0d0 * dlog(u1))
  theta = 2.0d0 * pi * u2
  z = r * dcos(theta)

  ! Scale and shift
  rnumber = mean + std * z

end subroutine gaussian_random


subroutine Initialize_RhoROCK
  implicit none

  real*8 :: rann

  do ic = 1, Nc

    call random_number(rann)
    !Rho(ic) = rann
    !Rho(ic) = rann*1.0d-1
    Rho(ic) = rann*beta

    call random_number(rann)
    !ROCK(ic) = rann
    !ROCK(ic) = rann*1.0d-1
    ROCK(ic) = rann*beta

    call random_number(rann)
    !Myosin(ic) = rann
    !Myosin(ic) = rann*1.0d-1
    Myosin(ic) = rann*beta

  end do

!  Rho = 0.0d0
!  ROCK = 0.0d0
!  Myosin = 0.0d0

end subroutine Initialize_RhoROCK


subroutine Solve_RhoROCK_Euler
  implicit none

  real*8 :: Rho_RHS(size(num)), ROCK_RHS(size(num)), Myosin_RHS(size(num))
  real*8 :: area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin
  real*8 :: Myosin_noise(Nc)


  Myosin_noise = 0.0d0
  if(if_myosin_noise)then
    call random_number(Myosin_noise)
  end if

  
  do ic = 1, Nc

    nn = num(ic)
    vx = v(1,inn(1:nn,ic))
    vy = v(2,inn(1:nn,ic))
    call CalculateArea(vx,vy,nn,area)
    area = abs(area)
    area_diff = area - Ao

    if (area_diff .lt. 0)then
      f_Rho = 0.0d0
    else
      f_Rho = A_Rho * area_diff**nhill /  (K_hill**nhill + area_diff**nhill) 
    end if


    f_ROCK = A_ROCK * Rho(ic)
    f_Myosin = A_Myosin * ROCK(ic)

    
    Rho_RHS(ic) = f_Rho * (1 - Rho(ic)) - D_Rho * Rho(ic)
    ROCK_RHS(ic) = f_ROCK * (1 - ROCK(ic)) - D_ROCK * ROCK(ic)
    Myosin_RHS(ic) = f_Myosin * (1 - Myosin(ic)) - D_Myosin * Myosin(ic)

  end do


  Rho(1:Nc) = Rho(1:Nc) +  dt * Rho_RHS(1:Nc)
  ROCK(1:Nc) = ROCK(1:Nc) + dt * ROCK_RHS(1:Nc)
  Myosin(1:Nc) = Myosin(1:Nc) + dt * Myosin_RHS(1:Nc) + & 
    sqrt(dt) * Myosin_noise_strength * Myosin_noise(1:Nc)



end subroutine Solve_RhoROCK_Euler

subroutine Solve_RhoROCK_RK4
  implicit none

  integer :: ic
  real*8 :: area, area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin

  ! RK4 intermediate arrays
  real*8 :: k1_Rho(Nc), k2_Rho(Nc), k3_Rho(Nc), k4_Rho(Nc)
  real*8 :: k1_ROCK(Nc), k2_ROCK(Nc), k3_ROCK(Nc), k4_ROCK(Nc)
  real*8 :: k1_M(Nc), k2_M(Nc), k3_M(Nc), k4_M(Nc)

  real*8 :: Rho_tmp(Nc), ROCK_tmp(Nc), M_tmp(Nc)

  real*8 :: Myosin_noise1(Nc), Myosin_noise2(Nc), propen_Myosin1(Nc), propen_Myosin2(Nc)
  real*8, parameter :: eps = 1.0d-3


  Myosin_noise1 = 0.0d0
  propen_Myosin1 = 0.0d0
  propen_Myosin2 = 0.0d0

  if(if_myosin_noise)then
    
    if(if_gaussian_noise)then
      do ic = 1, Nc
        call gaussian_random(Myosin_noise1(ic), 0.0d0, 1.0d0)
        call gaussian_random(Myosin_noise2(ic), 0.0d0, 1.0d0)
      end do
    else
      call random_number(Myosin_noise1)
      call random_number(Myosin_noise2)
      Myosin_noise1 = Myosin_noise1 - 0.5d0
      Myosin_noise2 = Myosin_noise2 - 0.5d0
    end if

    propen_Myosin1(1:Nc) = A_Myosin * ROCK(1:Nc) * (1.0d0 - Myosin(1:Nc))
    propen_Myosin2(1:Nc) = D_Myosin * Myosin(1:Nc)
  end if


  ! ============
  ! RK4 steps
  ! ============

  ! k1
  call compute_RHS(Rho, ROCK, Myosin, k1_Rho, k1_ROCK, k1_M)

  ! k2
  Rho_tmp  = Rho(1:Nc)   + 0.5d0 * dt * k1_Rho
  ROCK_tmp = ROCK(1:Nc)  + 0.5d0 * dt * k1_ROCK
  M_tmp    = Myosin(1:Nc) + 0.5d0 * dt * k1_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k2_Rho, k2_ROCK, k2_M)

  ! k3
  Rho_tmp  = Rho(1:Nc)   + 0.5d0 * dt * k2_Rho
  ROCK_tmp = ROCK(1:Nc)  + 0.5d0 * dt * k2_ROCK
  M_tmp    = Myosin(1:Nc) + 0.5d0 * dt * k2_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k3_Rho, k3_ROCK, k3_M)

  ! k4
  Rho_tmp  = Rho(1:Nc)   + dt * k3_Rho
  ROCK_tmp = ROCK(1:Nc)  + dt * k3_ROCK
  M_tmp    = Myosin(1:Nc) + dt * k3_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k4_Rho, k4_ROCK, k4_M)

  ! Final RK4 update
  Rho(1:Nc) = Rho(1:Nc)   + (dt/6.0d0) * (k1_Rho + 2*k2_Rho + 2*k3_Rho + k4_Rho) 
  ROCK(1:Nc) = ROCK(1:Nc)  + (dt/6.0d0) * (k1_ROCK + 2*k2_ROCK + 2*k3_ROCK + k4_ROCK)
  Myosin(1:Nc) = Myosin(1:Nc) + (dt/6.0d0) * (k1_M + 2*k2_M + 2*k3_M + k4_M) &
    + sqrt(dt) * Myosin(1:Nc) * (1.0d0 - Myosin(1:Nc)) * & 
    Myosin_noise_strength * Myosin_noise1(1:Nc)


!! -----------------------------------------------------------------------------------------
!  Myosin(1:Nc) = Myosin(1:Nc) + (dt/6.0d0) * (k1_M + 2*k2_M + 2*k3_M + k4_M) &
!    + sqrt(dt) * sqrt(Propen_Myosin1(1:Nc)) * Myosin_noise_strength * Myosin_noise1(1:Nc) &
!    - sqrt(dt) * sqrt(Propen_Myosin2(1:Nc)) * Myosin_noise_strength * Myosin_noise2(1:Nc)

!do ic = 1, Nc
!  if(Myosin(ic) .lt. eps) then
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!     + sqrt(dt) * sqrt(Propen_Myosin1(ic)) * Myosin_noise_strength * Myosin_noise1(ic) 

! elseif(Myosin(ic) .gt. 1.0d0 - eps)then
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!      - sqrt(dt) * sqrt(Propen_Myosin2(ic)) * Myosin_noise_strength * Myosin_noise2(ic)

!  else
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!     + sqrt(dt) * sqrt(Propen_Myosin1(ic)) * Myosin_noise_strength * Myosin_noise1(ic) &
!      - sqrt(dt) * sqrt(Propen_Myosin2(ic)) * Myosin_noise_strength * Myosin_noise2(ic)
!  end if

!  if (Myosin(ic) .lt. 0.0d0 .or. Myosin(ic) .gt. 1.0d0) then
!    print*, 'out of bounds', ic, Myosin(ic)
!  end if

!end do

!!--- Applying reflecting BC. 

!  do ic = 1, Nc
!    if(Myosin(ic) .gt. 1.0d0) then
!      Myosin(ic) = 2.0d0 - Myosin(ic)
!    end if
!    if(Myosin(ic) .lt. 0.0d0) then
!      Myosin(ic) =  - Myosin(ic)
!    end if
!  end do

!! -----------------------------------------------------------------------------------------

end subroutine Solve_RhoROCK_RK4

subroutine compute_RHS(Rho_in, ROCK_in, M_in, Rho_out, ROCK_out, M_out)
  ! BUGFIX (log.txt): missing implicit none let f_Rho/f_ROCK/f_Myosin fall
  ! back to Fortran's default implicit typing (single-precision real) instead
  ! of real*8, silently degrading precision in the RK4 solver path relative
  ! to the equivalent Euler path (Solve_RhoROCK_Euler), which declares them
  ! explicitly as real*8.
  implicit none
  real*8, intent(in)  :: Rho_in(Nc), ROCK_in(Nc), M_in(Nc)
  real*8, intent(out) :: Rho_out(Nc), ROCK_out(Nc), M_out(Nc)
  real*8 :: area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin

  do ic = 1, Nc
    nn = num(ic)
    vx = v(1, inn(1:nn,ic))
    vy = v(2, inn(1:nn,ic))

    call CalculateArea(vx,vy,nn,area)
    area = abs(area)
    area_diff = area - Ao

    ! Rho production
    if (area_diff .lt. 0.0d0) then
        f_Rho = 0.0d0
    else
        f_Rho = A_Rho * area_diff**nhill / (K_hill**nhill + area_diff**nhill)
    end if

    ! ROCK and Myosin production
    f_ROCK   = A_ROCK   * Rho_in(ic)
    f_Myosin = A_Myosin * ROCK_in(ic)

    ! RHS
    Rho_out(ic)  = f_Rho   * (1.0d0 - Rho_in(ic))  - D_Rho   * Rho_in(ic)
    ROCK_out(ic) = f_ROCK  * (1.0d0 - ROCK_in(ic)) - D_ROCK  * ROCK_in(ic)
    M_out(ic)    = f_Myosin * (1.0d0 - M_in(ic))    - D_Myosin * M_in(ic)
  end do

end subroutine compute_RHS

subroutine Solve_Myosin_Act_Contr

  implicit none
  integer :: ic
  real*8 :: rann

  do ic = 1,Nc 

    call random_number(rann)
    rann = rann - 0.5d0

    Myosin(ic) = Myosin(ic) - dt * ( Myosin(ic) - beta_0)/tau_contr &
      + sqrt(dt) * active_contr_strength * rann

  end do


end subroutine Solve_Myosin_Act_Contr
  




end module Force
