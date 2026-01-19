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

  Nc = Lx*Ly


  if(if_motility_gradient)then
    print*, 'Motility Gradient'
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

   if(if_active_contractility)then
     print*, 'Active Contractility ; ', 'strength:', active_contr_strength
   end if

   if(if_ABP)then
     print*, 'Active Brownian ; ', 'vo', vo, 'Dr', Dr
   end if


   if(if_RhoROCK)then
     print*, 'RhoROCK'
     call  Initialize_RhoROCK
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
      do iki = 1,100

        if(if_Do_T1)then
          call Do_T1
        end if

        if(if_Do_T2)then
          call Do_T2
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



    call Force_Calculation
    call Motile_Force_Calculation

    if(if_Fixed_boundary)then
      call Apply_Fixed_Boundary
    end if

     if(if_bottom_borders_fixed)then
       call Apply_bottom_border_Fixed
     end if

     if(if_top_borders_fixed)then
       call Apply_top_border_Fixed
     end if

     if(if_ABP)then
       call ABP_Force_Calculation
     end if




    v(1,:) = v(1,:) + dt * fxx(:)/eta + sqrt(dt) * fxx_ran(:)/eta + &
      sqrt(dt) * fxx_active_contr(:) / eta + &
      dt * fxx_ABP(:) / eta
      

    v(2,:) = v(2,:) + dt * fyy(:)/eta + sqrt(dt) * fyy_ran(:)/eta + &
      sqrt(dt) * fyy_active_contr(:) / eta + &
      dt * fyy_ABP(:) / eta

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


    write(*,'(A,I0,A,I0)', advance='no') achar(13)//'Step: ', it, '/', totT

!   write(*,*)'it= ', it
!   write(*,*)'tttttttttttttttttttttttttttttttt'


  end do ! ends it



!  write(*,*)"Total T1 count", Total_T1_count
!  write(*,*)"Total T2 count", Total_T2_count

  write(*,*)"**********End******************"

end program vertexmain
