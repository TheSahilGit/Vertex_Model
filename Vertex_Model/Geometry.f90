module Geometry 

  use array_info



  contains

   subroutine CalculateDistance(x1,y1,x2,y2,distance)
    implicit none
    real*8, intent(in) :: x1,y1,x2,y2
    real*8, intent(out) :: distance
    real*8 :: di


    di = (x1-x2)**2 + (y1-y2)**2
    distance = sqrt(di)

   end subroutine CalculateDistance

   subroutine CalculatePerimeter(x,y,dimn,perimeter)
       implicit none
       integer*4 :: i,j
       integer*4 :: dimn
       real*8, dimension(:), allocatable :: x,y
       real*8 :: perimeter
       real*8 :: d
   
   
   
   
       perimeter = 0.0d0
       d = 0.0d0
   
   
       do i=1,dimn
           j = i+1
           if(i==dimn)then
               j=1
           end if
           d = (x(i)-x(j))**2 + (y(i)-y(j))**2
           perimeter = perimeter + sqrt(d)
   
       end do

    end subroutine CalculatePerimeter

    subroutine CalculateArea(x,y,dimn,area)
      implicit none
      integer*4 :: i,j
      integer*4, intent(in) :: dimn
      real*8, dimension(:), allocatable :: x,y
      real*8 :: area
      real*8 :: sum1, sum2
      integer*4 :: n1,n2
  
  
  
  
  
      area = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
  
  
      do i=1,dimn
          j = i+1
          if(i==dimn)then
              j=1
          end if
  
          sum1 = sum1 + x(i)*y(j)
          sum2 = sum2 + x(j)*y(i)
  
          area = (sum1-sum2)/2.0
      end do

    end subroutine CalculateArea

    
    subroutine PerpendicularBisector(x,y,d,xnew, ynew)
      
      real*8, dimension(:), allocatable :: x,y, xnew, ynew
      real*8 :: xm, ym
      real*8 :: d !The ratio of length and the minimum length. 

      xm = (x(1) + x(2))/2.0
      ym = (y(1) + y(2))/2.0

      xnew(1) = -d * y(1) + d * ym + xm
      xnew(2) = -d * y(2) + d * ym + xm

      ynew(1) = d * x(1) - d * xm + ym
      ynew(2) = d * x(2) - d * xm + ym


    end subroutine PerpendicularBisector

      subroutine PerpendicularBisector2(x, y, d, xnew, ynew)
    implicit none

    ! Arguments
    real*8, intent(in) :: x(2), y(2)  ! Input coordinates
    real*8, intent(in) :: d           ! Scaling factor
    real*8, intent(out) :: xnew(2), ynew(2)  ! Output coordinates

    ! Local variables
    real*8 :: xm, ym

    ! Compute midpoints
    xm = (x(1) + x(2)) / 2.0d0
    ym = (y(1) + y(2)) / 2.0d0

    ! Calculate new coordinates based on perpendicular bisector
    xnew(1) = -d * y(1) + d * ym + xm
    xnew(2) = -d * y(2) + d * ym + xm

    ynew(1) = d * x(1) - d * xm + ym
    ynew(2) = d * x(2) - d * xm + ym

  end subroutine PerpendicularBisector2

   subroutine Get_Boundary_info

     use allocation

     implicit none
  
  
       ii = 0
       do i=1,Lx*Ly
         ii = ii+1
         mainarea(ii) = i
       end do
       ii = 0
       do i=1,Ly
         ii = ii+1
         leftP(ii) = i
       end do
       ii = 0
       do i=Lx*Ly-Ly+1,Lx*Ly
         ii = ii+1
         rightP(ii) = i
       end do
       ii = 0
       do i=Ly,Lx*Ly,Ly ! start,end,step
         ii = ii+1
         topP(ii) = i
       end do
       ii = 0
       do i=1,Lx*Ly-Ly+1,Ly
         ii = ii+1
         bottomP(ii) = i
       end do
       corners = [int(1),Ly,Lx*Ly-Ly+int(1),Lx*Ly]
  
       jj = 0
       do j=Lx*Ly+1,Lx*Ly+Ly
         jj = jj+1
         GrightP(jj) = j
       end do
       jj = 0
       do j = 3*Lx*Ly-Ly+1,3*Lx*Ly
         jj = jj+1
         GleftP(jj) = j
       end do
       jj = 0
       do j = 3*Lx*Ly+1,4*Lx*Ly-Ly+1,Ly
         jj = jj+1
         GtopP(jj) = j
       end do
       jj = 0
       do j = 4*Lx*Ly+Ly,5*Lx*Ly,Ly
         jj = jj+1
         GbottomP(jj) = j
       end do
       Gcorners = [5*Lx*Ly+int(1), 7*Lx*Ly+Ly, 9*Lx*Ly-Ly+int(1), 7*Lx*Ly]
  
       workingzone = [mainarea, GleftP, GrightP, GtopP, GbottomP, Gcorners]
  
       boundary = [leftP, rightP, topP, bottomP, corners]

        ii = 0
        ! Loop for inside1: Exclude 1 layer of boundary
        do i = 2, Lx - 1                ! Rows excluding first and last
          do j = 2, Ly - 1              ! Columns excluding first and last
            ii = ii + 1
            inside1(ii) = (i - 1) * Ly + j
          end do
        end do
 
        ii = 0
        ! Loop for inside2: Exclude 2 layers of boundary
        do i = 3, Lx - 2                ! Rows excluding first two and last two
          do j = 3, Ly - 2              ! Columns excluding first two and last two
            ii = ii + 1
            inside2(ii) = (i - 1) * Ly + j
          end do
        end do




 end subroutine Get_Boundary_info



 subroutine Find_boundary_dynamic

   use allocation
   implicit none

   logical :: is_boundary
   integer :: im
   integer:: boundary_temp(Lx*Ly*6)

   

  ! bottom_border, top_border, left_border, right_border

   real*8 :: maxVx, maxVy, minVx, minVy

     vertex_occurance_count = 0 

     do ic = 1,Lx*Ly
       nn = num(ic)
       do jc = 1,nn
         vertex_occurance_count(inn(jc,ic)) = vertex_occurance_count(inn(jc,ic)) + 1
       end do
     end do



     all_border_count = 0
     do jc = 1, size(vertex_occurance_count)

       if(vertex_occurance_count(jc).lt.3.and. &
         vertex_occurance_count(jc).gt.0)then

         all_border_count = all_border_count + 1

         all_borders(all_border_count) = jc

       end if
     end do

     minVx = minval(v(1,all_borders(1:all_border_count)))
     maxVx = maxval(v(1,all_borders(1:all_border_count)))
     minVy = minval(v(2,all_borders(1:all_border_count)))
     maxVy = maxval(v(2,all_borders(1:all_border_count)))


     bottom_border_count = 0
     top_border_count = 0
     left_border_count = 0
     right_border_count = 0

     do jc = 1, all_border_count

       if(v(2, all_borders(jc)).lt.minVy+1.0d0.and.v(2, all_borders(jc)).gt.minVy-1.0d0)then
         bottom_border_count = bottom_border_count + 1
         bottom_border(bottom_border_count) = all_borders(jc)
       end if

       if(v(2, all_borders(jc)).lt.maxVy+1.0d0.and.v(2, all_borders(jc)).gt.maxVy-1.0d0)then
         top_border_count = top_border_count + 1
         top_border(top_border_count) = all_borders(jc)
       end if

       if(v(1, all_borders(jc)).lt.minVx+1.0d0.and.v(1, all_borders(jc)).gt.minVx-1.0d0)then
         left_border_count = left_border_count + 1
         left_border(left_border_count) = all_borders(jc)
       end if

       if(v(1, all_borders(jc)).lt.maxVx+1.0d0.and.v(1, all_borders(jc)).gt.maxVx-1.0d0)then
         right_border_count = right_border_count + 1
         right_border(right_border_count) = all_borders(jc)
       end if
     end do


      boundary_count = 0
      do ic = 1, Lx * Ly
        nn = num(ic)
        is_boundary = .false.

        do jc = 1, nn
          do im = 1, all_border_count

          if (all_borders(im) == inn(jc,ic)) then
            is_boundary = .true.
            boundary_count = boundary_count + 1
            boundary_temp(boundary_count) = ic
            exit  ! If found, mark the cell as a boundary
          end if

        end do
        end do


      end do

     call find_unique_sorted(boundary_temp, boundary_count, boundary, boundary_occ)




 end subroutine Find_boundary_dynamic




end module Geometry
