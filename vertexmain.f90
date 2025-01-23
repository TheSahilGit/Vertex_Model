! Here what I am not using any fortran library or subroutine to create the initial voronoi structure. 
! As, I need voronoi only to create the initial sturcture, I do it using matlab ( It is easier their) and write those
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


  call read_input
  call allocate_arrays
  call read_data   ! Initialization


  Total_T1_count = 0
  Total_T2_count = 0


  if(if_motility_gradient)then
     call Give_Motility_Gradient
  else
     mot = etas_max
  end if
  mot0 = mot


  do it = 1,totT

    
    call Calculate_Energy

!*********************************************************************
    if(modulo(it,it_dump).eq.0.or.it.eq.1.or.it.eq.2)then
      call write_output
    end if

!******************************************************************
   



    if(modulo(it,T1_time_interval).eq.0.0d0.or.it.eq.1)then
!    if(it==1)then
      do iki = 1,100
        call Do_T1
      end do
    end if


    if(modulo(it,T2_time_interval).eq.0.0d0)then
      do iki= 1, 1
        call Do_T2
      end do
    end if





    call Force_Calculation
    call Motile_Force_Calculation

    if(if_Fixed_boundary)then
      call Apply_Fixed_Boundary
    end if

     if(if_bottom_borders_fixed)then
       call Apply_bottom_border_Fixed
     end if




    v(1,:) = v(1,:) + dt * fxx(:)/eta + sqrt(dt) * fxx_ran(:)/eta
    v(2,:) = v(2,:) + dt * fyy(:)/eta + sqrt(dt) * fyy_ran(:)/eta

    
    if(if_Perturb_tissue)then
      call Apply_perturbation
    end if

    if(if_Shear_tissue)then
      call ShearTissue
      call Calculate_StressTensor
    end if


   write(*,*)'it= ', it
!   write(*,*)'tttttttttttttttttttttttttttttttt'


  end do ! ends it



  write(*,*)"Total T1 count", Total_T1_count
  write(*,*)"Total T2 count", Total_T2_count

  write(*,*)"**********End******************"

end program vertexmain
