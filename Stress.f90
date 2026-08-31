module Stress

  use allocation 
  use System_Info


  contains

   subroutine Calculate_StressTensor
 
 
       implicit none
       integer :: cellNo
       real*8 :: term1, term2, dx, dy, len_d
 
     !  call Get_Boundary_info
       
       radius_from_core = sqrt(dble(Lx*Lx) + dble(Ly*Ly))/5.0d0
       !radius_from_core = 4.0d0 * sqrt(Ao)

       call Get_Cells_Within_Radius

!       open(1211, file='inside.dat', status='unknown')
!         write(1211,*)n_inside, inside_cells(1:n_inside)
!       close(1211)
 
       ! BUGFIX (log.txt): beta/gamm nondimensionalization moved to read_input
       ! (allocation.f90), done once instead of every call to this subroutine
       ! -- see log.txt for why repeating it here compounded beta/gamm.

       totalarea = 0.0d0
       TotalSigma = 0.0d0

 
       do ic = 1,n_inside
          cellNo = inside_cells(ic)
          nn = num(cellNo)
 
          vx = v(1,inn(1:nn,cellNo))
          vy = v(2,inn(1:nn,cellNo))
 
          call CalculateArea(vx,vy,nn,area)
          area = abs(area)
          totalarea = totalarea + area
       end do
 
 
       do ic = 1, n_inside
         cellNo = inside_cells(ic)
         nn = num(cellNo)
 
         vx = v(1,inn(1:nn,cellNo))
         vy = v(2,inn(1:nn,cellNo))
 
         call CalculateArea(vx,vy,nn,area)
         area = abs(area)
         call CalculatePerimeter(vx,vy,nn,perimeter)
 
         term1 =  2.0d0 * lambda * (area - Ao)
         term2 = (2.0d0 * beta* perimeter + gamm)/(2.0d0*area)

         ! BUGFIX (log.txt): sigma is a global (allocation.f90) that was never
         ! reset per cell, so each cell's edge-sum accumulated on top of the
         ! *previous cell's already fully-scaled* tensor -- corrupting
         ! TotalSigma/ShearStress for every cell after the first, every call.
         sigma = 0.0d0

         do jc = 1,nn
           jp = jc+1
           if (jc.eq.nn)then
             jp = 1
           end if
           dx = vx(jp) - vx(jc)
           dy = vy(jp) - vy(jc)
           len_d = sqrt(dx*dx + dy*dy)
 
           sigma(1,1) = sigma(1,1) + dx*dx/len_d
           sigma(1,2) = sigma(1,2) + dx*dy/len_d
           sigma(2,1) = sigma(2,1) + dy*dx/len_d
           sigma(2,2) = sigma(2,2) + dy*dy/len_d
 
 
         end do
 
           sigma(1,1) = term1  + term2 * sigma(1,1)
           sigma(1,2) = term2 * sigma(1,2)
           sigma(2,1) = term2 * sigma(2,1)
           sigma(2,2) = term1 + term2 * sigma(2,2)
 
 
 
 
        TotalSigma = TotalSigma + sigma * area/totalarea
 
 
        end do
 
 
        ShearStress(it) = TotalSigma(1,2)
 
    end subroutine Calculate_StressTensor



    subroutine ShearTissue

    implicit none

    real*8 :: comb

     ! write(*,*)if_Shearing, shearStrength

     if(if_Sudden_Shearing)then
       if(it.eq.sudden_shearWhen)then
           
          v(1,:) = v(1,:) + sudden_shearStrength * v(2,:)
          v(2,:) = v(2,:)

          if(if_bottom_borders_fixed)then
            call Find_boundary_dynamic
            v(1, bottom_border(1:bottom_border_count)) = &
              v(1,bottom_border(1:bottom_border_count)) - &
              sudden_shearStrength * v(2,bottom_border(1:bottom_border_count))
          end if

       end if
     end if 

     if(if_Oscillatory_Shearing)then
       if(it.gt.Oscl_shearWhen)then


         strainRate = Oscl_shearStrength * Oscl_freq_wo* cos(Oscl_freq_wo * it * dt)

         v(1,:) = v(1,:) + dt * strainRate * v(2,:)      
    
       end if

     end if




    end subroutine ShearTissue





end module Stress
