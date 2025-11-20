module System_Info

  use allocation
  use Geometry

  contains
    subroutine Calculate_Energy
      
      real*8 :: En
      real*8 :: En_Area, En_Peri, En_Line




      En = 0.0d0

      do ic = 1, Nc !Lx*Ly

        nn = num(ic)
        
        vx = v(1, inn(1:nn,ic))
        vy = v(2, inn(1:nn,ic))

        call CalculateArea(vx, vy, nn, area)
        call CalculatePerimeter(vx, vy, nn, perimeter)
        area = abs(area)


        En_Area =  lambda * (area - Ao)**2
        En_Peri =  beta * (perimeter - Co)**2
        En_Line =  gamm * perimeter
        
        En = En + En_Area + En_Peri + En_Line
     
       end do

       Energy(it) = En


    end subroutine Calculate_Energy


   subroutine CellCentre
!     use Geometry
!     use allocation

     do ic = 1, Nc !Lx*Ly

       vx = v(1, inn(1:num(ic),ic))
       vy = v(2, inn(1:num(ic),ic))

       cellcen(ic,1) = sum(vx)/dble(size(vx))
       cellcen(ic,2) = sum(vy)/dble(size(vy))

      end do

      global_cellCenX = sum(cellcen(:,1))/dble(Nc) ! Lx*Ly
      global_cellCenY = sum(cellcen(:,2))/dble(Nc) ! Lx*Ly

      



   end subroutine CellCentre




end module System_Info
