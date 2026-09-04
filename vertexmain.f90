! Here what I am not using any fortran library or subroutine to create the initial voronoi structure. 
! As, I need voronoi only to create the initial sturcture, I do it using matlab ( It is easier there) and write those
! data in a '*.dat' file, and then load those data in this fortran code. 
! It would not slow down the speed as we are gonna do it only once at the start and use the data file for the rest.
! (The time taken to create voronoi diagram with 100x100 points is approx 30 sec.) 
! REMEMBER : If we want to make changes like perfect honeycomb to randomly structured points, those we'll have to 
! do in the matlab code 'Main.m' . 



program vertexmain
  use allocation
  use array_info
  use System_Info
  use T1_Swap
  use T2_Swap
  use Force
  use Stress
  use Proliferation


  call read_input
  call allocate_arrays
  call read_data   ! Initialization



  Total_T1_count = 0
  Total_T2_count = 0

  ! BUGFIX (log.txt): Nc is now derived from the loaded num() array at the
  ! end of read_data (allocation.f90), correctly for both nrun=1 and
  ! nrun=2 -- see that comment for why hardcoding it to Lx*Ly here broke
  ! nrun=2 restarts taken after any division/T2 event.

  if(if_motility)then
    print*, "Apolar Cell Motility, ", ' etas_max', etas_max
  end if

  if(if_motility_gradient)then
    if(if_motility_Eulerian)then
      print*, 'Motility Gradient (Eulerian -- recomputed from current position every step)'
    else
      print*, 'Motility Gradient (Lagrangian -- assigned once, carried per vertex)'
    end if
    call Give_Motility_Gradient
  elseif(if_motility_hotspot)then
    print*, 'Motility Hotspot'
    call Give_Motility_Hotspot
  else
     mot = etas_max
  end if
  mot0 = mot


  if(if_Fixed_boundary)then
    print*, "Fixed boundary"
  end if

   if(if_bottom_borders_fixed)then
     print*, "Bottom border fixed"
   end if

   if(if_top_borders_fixed)then
     print*, "Top border fixed"
   end if

   if(if_cell_division)then
     print*, "Cell division true. A_cut", area_0
   end if

   if(if_active_contractility)then
     print*, 'Active Contractility ; ', 'tau', tau_contr, 'strength:', active_contr_strength
     Myosin = beta_0
   end if

   if(if_ABP)then
     print*, 'Active Brownian ; ', 'vo', vo, 'Dr', Dr
   end if


   if(if_RhoROCK)then
     print*, 'RhoROCK'
     call  Initialize_RhoROCK
   end if

  if(if_polar_motility)then
    print*, 'Polar Motility ; ' , 'Strength', polar_motility_strength
  end if




  do it = 1,totT

    
    call Calculate_Energy

!*********************************************************************
    if(modulo(it,it_dump).eq.0.or.it.eq.1.or.it.eq.2)then
      call write_output
    end if

!******************************************************************
   

   if(modulo(it,T2_time_interval).eq.0.0d0)then
     do iki= 1, 1
!       call Do_T2
     end do
   end if


    if(modulo(it,T1_time_interval).eq.0.0d0.or.it.eq.1)then
!    if(it==1)then
      do iki = 1,500

        if(if_Do_T1)then
          call Do_T1
        end if

        if(if_Do_T2)then
          call Do_T2
        end if

        ! OPTIMIZATION (log.txt): find_T1/find_T2 (called at the top of
        ! Do_T1/Do_T2) do a full O(Nc) mesh scan every one of these 500
        ! iterations, even after no short edges/small triangles remain --
        ! count_T1/count_T2 are freshly recomputed by that scan on every
        ! call. Once both are simultaneously zero, nothing in the mesh can
        ! change for the rest of this loop (no candidates -> no possible
        ! mutation -> every remaining full rescan would deterministically
        ! find the same empty result), so stopping here skips zero
        ! legitimate T1/T2 events -- only the wasted rescans that would
        ! have found nothing anyway.
        if((.not.if_Do_T1 .or. count_T1==0) .and. (.not.if_Do_T2 .or. count_T2==0)) then
          exit
        end if

      end do
    end if





    if(if_RhoROCK)then
      if(if_RK4)then
        call Solve_RhoROCK_RK4
      else
        call Solve_RhoROCK_Euler
      end if
    end if

    if(if_active_contractility)then
      call Solve_Myosin_Act_Contr
    end if



    call Force_Calculation

    if(if_motility)then
      ! Eulerian motility (log.txt, user request): if_motility_gradient's
      ! initial assignment (above, before this loop) is a one-time
      ! Lagrangian tag -- correct per vertex ID forever after, but blind to
      ! where that vertex has since moved. if_motility_Eulerian re-derives
      ! mot from every vertex's CURRENT position each step instead, so the
      ! gradient stays anchored to physical location rather than drifting
      ! with the tissue. Cheap enough to call unconditionally every step
      ! after the Give_Motility_Gradient optimization (Force.f90) -- see
      ! that subroutine's comment. NOTE: if_motility_decay's effect is
      ! largely moot combined with this -- decay modifies mot in place
      ! after this call, but the very next step's recompute discards it,
      ! so decay never accumulates under Eulerian mode.
      if(if_motility_gradient .and. if_motility_Eulerian)then
        call Give_Motility_Gradient
      end if
      call Motile_Force_Calculation
    end if

     if(if_ABP)then
       call ABP_Force_Calculation
     end if

    if(if_polar_motility)then
      call Polar_Motile_Force_Calculation
    end if

    !! -- Apply Boundary conditions --

    if(if_Fixed_boundary)then
      call Apply_Fixed_Boundary
    end if

     if(if_bottom_borders_fixed)then
       call Apply_bottom_border_Fixed
     end if

     if(if_top_borders_fixed)then
       call Apply_top_border_Fixed
     end if
    
   !! ---------------------------------------



    v(1,:) = v(1,:) + dt * fxx(:)/eta + sqrt(dt) * fxx_ran(:)/eta + &
      dt * fxx_ABP(:) / eta + &
      sqrt(dt) * fxx_Polar(:)/eta
      

    v(2,:) = v(2,:) + dt * fyy(:)/eta + sqrt(dt) * fyy_ran(:)/eta + &
      dt * fyy_ABP(:) / eta + &
      sqrt(dt) * fyy_Polar(:)/eta

    if(if_ABP)then
      theta_ABP(:)  = theta_ABP(:) + dsqrt(dt) * rot_noise(:)
    end if


    
    if(if_Perturb_tissue)then
      call Apply_perturbation
    end if

    if(if_Shear_tissue)then
      call ShearTissue
      call Calculate_StressTensor
    end if



    if(if_limb_force)then
      call Apply_Limb_Force
    end if

    if(if_cell_division)then
      call Do_Proliferation
    end if

    if(if_squeeze_tissue)then
      if(it.eq.squeeze_when)then
        call Squeeze_Tissue
      end if
    end if

    ! PBC (log.txt): wrap every vertex back into the canonical
    ! [0,Lx_box) x [0,Ly_box) box once per step, after every
    ! position-modifying step above has run. T1/T2/division already wrap
    ! their own newly-written/moved vertex positions individually at the
    ! point of writing (Wrap_Position calls in T1_swap.f90/T2_swap.f90/
    ! Proliferation.f90) -- wrapping again here is idempotent/harmless for
    ! those, and catches the ordinary force-integration position update
    ! (and Apply_perturbation/ShearTissue, if ever combined with if_PBC)
    ! in one place. A no-op when if_PBC is false.
    if(if_PBC)then
      v(1,:) = modulo(v(1,:), Lx_box)
      v(2,:) = modulo(v(2,:), Ly_box)
    end if

    write(*,'(A,I0,A,I0)', advance='no') achar(13)//'Step: ', it, '/', totT

!   write(*,*)'it= ', it
!   write(*,*)'tttttttttttttttttttttttttttttttt'


  end do ! ends it



!  write(*,*)"Total T1 count", Total_T1_count
!  write(*,*)"Total T2 count", Total_T2_count

  write(*,*)"**********End******************"

end program vertexmain
