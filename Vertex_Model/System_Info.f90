module System_Info

  use allocation
  use Geometry

  contains
    subroutine Calculate_Energy
      
      real*8 :: En
      real*8 :: En_Area, En_Peri, En_Line




      En = 0.0d0

      do ic = 1, Lx*Ly

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
     use allocation

     do ic = 1,Lx*Ly

       vx = v(1, inn(1:num(ic),ic))
       vy = v(2, inn(1:num(ic),ic))

       cellcen(ic,1) = sum(vx)/dble(size(vx))
       cellcen(ic,2) = sum(vy)/dble(size(vy))

      end do

      global_cellCenX = sum(cellcen(:,1))/dble(Lx*Ly)
      global_cellCenY = sum(cellcen(:,2))/dble(Lx*Ly)

      



   end subroutine CellCentre


!   subroutine MeanSqDisp(Lx_in,Ly_in,v_in,inn_in,num_in,cellCentInit,avgdisp)
!!     use allocation
!      
!     integer*4 :: ic, jc, nnm
!      integer*4, dimension(:), allocatable :: num_in
!      real*8, dimension(:,:), allocatable:: v_in
!      integer*4, dimension(:,:), allocatable :: inn_in
!      integer*4 :: Lx_in,Ly_in, working_L
!      real*8, allocatable, dimension(:) :: cellcen, cellCentInit
!      real*8, allocatable, dimension(:) :: disp
!      real*8 :: avgdisp
!
!
!
!!      write(*,*)'cellcent = ,', cellcen
!      
!      working_L = Lx_in*Ly_in - 4*Lx_in - 4*Ly_in + 16
!      
!      allocate(cellcen(working_L))
!      allocate(disp(working_L))
!      
!      call CellCentre(Lx_in,Ly_in,v_in,inn_in,num_in,cellcen) 
!
!      do i = 1, working_L
!        disp(i) = (abs(cellcen(i)-cellCentInit(i)))**2
!      end do
!
!      avgdisp = sum(disp)/(working_L)
!
!
!
!   end subroutine MeanSqDisp


end module System_Info
