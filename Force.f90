module Force
  use allocation
  use Geometry

  contains


  subroutine Force_Calculation

    implicit none

    integer :: prev_idx, next_idx
    real*8 :: len_d, dx, dy
    real*8 :: grad_perimeter_X, grad_perimeter_Y, grad_area_X, grad_area_Y


    fxx = 0.0d0
    fyy = 0.0d0



    beta = beta/(lambda*Ao)
    gamm = gamm/(lambda*(Ao)**1.5)

    do ic = 1, Nc !Lx*Ly
      
      nn = num(ic)

      vx = v(1,inn(1:nn,ic))
      vy = v(2,inn(1:nn,ic))
      


      call CalculateArea(vx,vy,nn,area)
      call CalculatePerimeter(vx,vy,nn,perimeter)
      area = abs(area)
      
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


        dx = vx(jc) - vx(prev_idx)
        dy = vy(jc) - vy(prev_idx)
        len_d = dsqrt(dx**2 + dy**2) 

        grad_perimeter_X = dx/len_d
        grad_perimeter_Y = dy/len_d


        dx = vx(jc) - vx(next_idx)
        dy = vy(jc) - vy(next_idx)
        len_d = dsqrt(dx**2 + dy**2) 

        grad_perimeter_X = grad_perimeter_X + dx/len_d
        grad_perimeter_Y = grad_perimeter_Y + dy/len_d


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

    real*8 :: rann
  

    if(if_motility_decay)then
     ! mot = mot0*exp(-it*dt/motility_decay_timeScale)

     mot = mot - dt * mot/motility_decay_timeScale
    end if

    do ic = 1,size(mot)
      call random_number(rann)
      fxx_ran(ic) = sqrt(2*mot(ic)*eta)*(2.0d0 * rann-1)
      call random_number(rann)
      fyy_ran(ic) = sqrt(2*mot(ic)*eta)*(2.0d0 * rann-1)
    end do


  end subroutine Motile_Force_Calculation


  subroutine Give_Motility_Gradient

    implicit none


    integer :: ip, idx
    real(8) :: vy, lowpatch, highpatch, meanpatch
    real(8) :: patchwidth
    integer :: verryCount(size(v, 2)), verry(size(v, 2))
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

  integer, parameter :: number_of_hotspot = 4
  integer, dimension(number_of_hotspot) :: hotspot_location  ! Cell index
  real*8, dimension(number_of_hotspot) :: xCM, yCM

  real*8 :: sigma_hotspot
  real*8 :: xij, yij, rij


  hotspot_location = (/ 4, Ly-4, Lx*Ly-Ly+1+4,  Lx*Ly-4 /)

  do jp = 1, number_of_hotspot
    ip = hotspot_location(jp)
    xCM(jp) = sum(v(1,inn(1:num(ip),ip)))/dble(num(ip))
    yCM(jp) = sum(v(2,inn(1:num(ip),ip)))/dble(num(ip))
  end do

!  print*, hotspot_location
!  print*, xCM
!  print*, yCM

  sigma_hotspot = 3.0d0

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
   end do


end subroutine Apply_Fixed_Boundary

subroutine Apply_bottom_border_Fixed
  implicit none

!  call Find_boundary_dynamic

  fxx(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_ran(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_ran(bottom_border(1:bottom_border_count)) = 0.0d0


end subroutine  Apply_bottom_border_Fixed

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

      if(ic .lt. int(Lx*Ly/2))then   ! Cheap way
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (-1.0d0)
      else 
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (1.0d0)
      end if

    end do







end subroutine Apply_Limb_Force



end module Force
