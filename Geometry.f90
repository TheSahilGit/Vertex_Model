module Geometry 

  use array_info
  use allocation



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

     do ic = 1,Nc !Lx*Ly
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
      do ic = 1, Nc !Lx * Ly
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



   subroutine Find_Principle_Axis(vx, vy, nn, x1, y1, x2, y2)
    implicit none
    integer, intent(in) :: nn
    real(8), intent(in) :: vx(nn), vy(nn)
    real(8), intent(out) :: x1, y1, x2, y2
  
    ! local variables
    integer :: i
    real(8) :: xm, ym
    real(8) :: Ixx, Iyy, Ixy
    real(8) :: theta, length, cos_t, sin_t
    real(8) :: x_min, x_max, proj
  
    ! ---- Compute centroid ----
    xm = sum(vx) / nn
    ym = sum(vy) / nn
  
    ! ---- Compute second moments (covariance matrix elements) ----
    Ixx = 0.0d0
    Iyy = 0.0d0
    Ixy = 0.0d0
    do i = 1, nn
       Ixx = Ixx + (vx(i) - xm)**2
       Iyy = Iyy + (vy(i) - ym)**2
       Ixy = Ixy + (vx(i) - xm)*(vy(i) - ym)
    end do
  
    ! ---- Compute principal direction ----
    theta = 0.5d0 * atan2(2.0d0 * Ixy, (Ixx - Iyy))
    cos_t = cos(theta)
    sin_t = sin(theta)
  
    ! ---- Project all points onto principal axis to find extents ----
    x_min = 1.0d99
    x_max = -1.0d99
    do i = 1, nn
       proj = (vx(i) - xm)*cos_t + (vy(i) - ym)*sin_t
       x_min = min(x_min, proj)
       x_max = max(x_max, proj)
    end do
  
    ! ---- Endpoints of principal axis ----
    x1 = xm + x_min*cos_t
    y1 = ym + x_min*sin_t
    x2 = xm + x_max*cos_t
    y2 = ym + x_max*sin_t
  
  end subroutine Find_Principle_Axis

subroutine Find_Bisector_Intersections(vx, vy, nn, x1, y1, x2, y2, xi, yi, &
    n_intersections, inn_input, ic, idx_pair)
  implicit none
  integer, intent(in) :: nn, ic
  real(8), intent(in) :: vx(nn), vy(nn)
  real(8), intent(in) :: x1, y1, x2, y2
  integer, intent(in) :: inn_input(nn)  
  real(8), intent(out) :: xi(nn), yi(nn)
  integer, intent(out) :: n_intersections
  integer, intent(out) :: idx_pair(2,2)

  integer :: i
  real(8) :: mx, my           ! midpoint
  real(8) :: dx, dy           ! principal axis direction
  real(8) :: nx, ny           ! bisector direction (normal)
  real(8) :: x3, y3, x4, y4   ! edge endpoints
  real(8) :: denom, ua, numer
  real(8) :: xtemp, ytemp
  real(8), parameter :: eps = 1.0d-12

  n_intersections = 0

  ! midpoint of principal axis
  mx = 0.5d0*(x1 + x2)
  my = 0.5d0*(y1 + y2)

  ! direction of principal axis
  dx = x2 - x1
  dy = y2 - y1

  ! perpendicular (normal) => bisector direction
  nx = -dy
  ny =  dx

  ! If the axis endpoints are identical (degenerate), return 0
  if (abs(nx) < eps .and. abs(ny) < eps) return

  do i = 1, nn
     x3 = vx(i)
     y3 = vy(i)
     if (i < nn) then
        x4 = vx(i+1)
        y4 = vy(i+1)
     else
        x4 = vx(1)
        y4 = vy(1)
     end if

     ! denom = (edge) x (bisector)  (2D scalar cross)
     denom = (x4 - x3)*ny - (y4 - y3)*nx

     if (abs(denom) > eps) then
        ! NOTE THE FIXED SIGN: numerator uses (mx - x3, my - y3)
        numer = (mx - x3)*ny - (my - y3)*nx
        ua = numer / denom

        if (ua >= 0.d0 .and. ua <= 1.d0) then
           xtemp = x3 + ua*(x4 - x3)
           ytemp = y3 + ua*(y4 - y3)
           n_intersections = n_intersections + 1
           xi(n_intersections) = xtemp
           yi(n_intersections) = ytemp
           
           if (n_intersections <= 2) then
              idx_pair(n_intersections,1) = inn_input(i)
              if (i < nn) then
                 idx_pair(n_intersections,2) = inn_input(i+1)
              else
                 idx_pair(n_intersections,2) = inn_input(1)
              end if
           end if

        end if

     else
        ! denom ~ 0 => edge is parallel to bisector direction.
        ! Could be collinear (infinite intersections) if midpoint lies on edge line.
        ! We'll detect if midpoint projects onto the edge line and, if so, record endpoints
        ! (or skip). Here we check collinearity and, if midpoint lies between endpoints,
        ! we record that projection point.
        numer = (mx - x3)*(y4 - y3) - (my - y3)*(x4 - x3)
        if (abs(numer) <= eps) then
           ! midpoint lies on the same infinite line as edge -> treat carefully
           ! project midpoint onto the edge segment and accept if inside segment
           ! compute t so that P = P3 + t*(P4-P3)
           if (abs(x4 - x3) >= abs(y4 - y3)) then
              ! avoid division by tiny x4-x3
              ua = (mx - x3) / (x4 - x3)
           else
              ua = (my - y3) / (y4 - y3)
           end if

           if (ua >= 0.d0 .and. ua <= 1.d0) then
              xtemp = x3 + ua*(x4 - x3)
              ytemp = y3 + ua*(y4 - y3)
              n_intersections = n_intersections + 1
              xi(n_intersections) = xtemp
              yi(n_intersections) = ytemp
              
              if (n_intersections <= 2) then
                idx_pair(n_intersections,1) = inn_input(i)
                if (i < nn) then
                   idx_pair(n_intersections,2) = inn_input(i+1)
                else
                   idx_pair(n_intersections,2) = inn_input(1)
                end if
              end if

           end if
        end if
     end if
  end do

end subroutine Find_Bisector_Intersections



subroutine SplitPolygon(inn_ic, nn, pair1, pair2, new_idx1, new_idx2, &
    inn_new1, n1, inn_new2, n2)
  implicit none
  integer, intent(in)  :: nn
  integer, intent(in)  :: inn_ic(nn)
  integer, intent(in)  :: pair1(2), pair2(2)
  integer, intent(in)  :: new_idx1, new_idx2
  integer, intent(out) :: inn_new1(nn+2), inn_new2(nn+2)
  integer, intent(out) :: n1, n2

  integer :: i, k
  integer :: posA, posB, posC, posD
  logical, allocatable :: in_between(:)
  integer :: idx

  ! --- initialize ---
  n1 = 0
  n2 = 0
  inn_new1 = 0
  inn_new2 = 0

  ! --- find positions of the chord vertices in inn_ic ---
  posA = -1; posB = -1; posC = -1; posD = -1
  do i = 1, nn
     if (inn_ic(i) == pair1(1)) posA = i
     if (inn_ic(i) == pair1(2)) posB = i
     if (inn_ic(i) == pair2(1)) posC = i
     if (inn_ic(i) == pair2(2)) posD = i
  end do

  if (posA == -1 .or. posB == -1 .or. posC == -1 .or. posD == -1) then
     print *, 'SplitPolygon error: chord vertices not found in inn_ic'
     return
  end if

  ! allocate boolean mask for vertices strictly between pair1(1) and pair2(2)
  allocate(in_between(nn))
  in_between = .false.

  ! mark indices strictly between posA and posD when walking forward:
  ! start at next(posA) and stop when you reach posD (do NOT mark posD)
  k = posA
  do
     k = k + 1
     if (k > nn) k = 1
     if (k == posD) exit
     in_between(k) = .true.
  end do

  ! -------------------------
  ! Build inn_new1:
  ! - take all vertices in order 1..nn except those marked 'in_between'
  ! - insert new_idx1,new_idx2 immediately after pair1(1) (posA)
  ! This yields: [inn_ic(1..posA), new1, new2, inn_ic(posD..end)] (no duplication)
  ! -------------------------
  n1 = 0
  do i = 1, nn
     if (.not. in_between(i)) then
        n1 = n1 + 1
        if (n1 > nn+2) then
           print *, 'SplitPolygon: overflow constructing inn_new1'
           deallocate(in_between)
           return
        end if
        inn_new1(n1) = inn_ic(i)
        if (i == posA) then
           ! insert new vertices right after pair1(1)
           n1 = n1 + 1
           inn_new1(n1) = new_idx1
           n1 = n1 + 1
           inn_new1(n1) = new_idx2
        end if
     end if
  end do

  ! -------------------------
  ! Build inn_new2:
  ! - vertices from after pair1(2) up to pair2(1) (forward, inclusive of posC)
  ! - then new_idx2,new_idx1 (reversed), then pair1(2) to close
  ! -------------------------
  n2 = 0
  idx = posB + 1
  if (idx > nn) idx = 1
  do
     n2 = n2 + 1
     if (n2 > nn+2) then
        print *, 'SplitPolygon: overflow constructing inn_new2'
        deallocate(in_between)
        return
     end if
     inn_new2(n2) = inn_ic(idx)
     if (idx == posC) exit
     idx = idx + 1
     if (idx > nn) idx = 1
  end do

  ! append new indices in reversed order then pair1(2)
  n2 = n2 + 1
  inn_new2(n2) = new_idx2
  n2 = n2 + 1
  inn_new2(n2) = new_idx1
  n2 = n2 + 1
  inn_new2(n2) = pair1(2)

  deallocate(in_between)

  return



end subroutine SplitPolygon


subroutine UpdateNeighborPolygons(inn_old1, n1, inn_old2, n2, &
                                  edge1, edge2, new_idx1, new_idx2, &
                                  inn_new1, n1_new, inn_new2, n2_new)
  implicit none
  integer, intent(in) :: n1, n2, new_idx1, new_idx2
  integer, intent(in) :: inn_old1(:), inn_old2(:)
  integer, intent(in) :: edge1(2), edge2(2)
  integer, intent(out) :: inn_new1(size(inn_old1)+2), inn_new2(size(inn_old2)+2)
  integer, intent(out) :: n1_new, n2_new

  integer :: i, a, b
  integer :: tmp(size(inn_old1)+2), tmp2(size(inn_old2)+2)
  integer :: m, m2

  ! init outputs
  n1_new = 0
  n2_new = 0
  inn_new1 = 0
  inn_new2 = 0

  ! -------------------------------
  ! Update first neighbor: check BOTH edges
  ! -------------------------------
  m = 0
  do i = 1, n1
     ! copy current vertex
     m = m + 1
     tmp(m) = inn_old1(i)

     ! determine next vertex (cyclic)
     if (i == n1) then
        a = inn_old1(i)
        b = inn_old1(1)
     else
        a = inn_old1(i)
        b = inn_old1(i+1)
     end if

     ! if this edge matches edge1 (either orientation) insert new_idx1 AFTER a
     if ( (a == edge1(1) .and. b == edge1(2)) .or. (a == edge1(2) .and. b == edge1(1)) ) then
        m = m + 1
        tmp(m) = new_idx1
     end if

     ! if this edge matches edge2 (either orientation) insert new_idx2 AFTER a
     if ( (a == edge2(1) .and. b == edge2(2)) .or. (a == edge2(2) .and. b == edge2(1)) ) then
        m = m + 1
        tmp(m) = new_idx2
     end if
  end do

  n1_new = m
  inn_new1(1:n1_new) = tmp(1:n1_new)

  ! -------------------------------
  ! Update second neighbor: check BOTH edges
  ! -------------------------------
  m2 = 0
  do i = 1, n2
     m2 = m2 + 1
     tmp2(m2) = inn_old2(i)

     if (i == n2) then
        a = inn_old2(i)
        b = inn_old2(1)
     else
        a = inn_old2(i)
        b = inn_old2(i+1)
     end if

     if ( (a == edge1(1) .and. b == edge1(2)) .or. (a == edge1(2) .and. b == edge1(1)) ) then
        m2 = m2 + 1
        tmp2(m2) = new_idx1
     end if

     if ( (a == edge2(1) .and. b == edge2(2)) .or. (a == edge2(2) .and. b == edge2(1)) ) then
        m2 = m2 + 1
        tmp2(m2) = new_idx2
     end if
  end do

  n2_new = m2
  inn_new2(1:n2_new) = tmp2(1:n2_new)

end subroutine UpdateNeighborPolygons





subroutine ArrangeVertices(v, v_dim1, v_dim2, inn_new1, n1, inn_new2, n2)
  implicit none

  ! Inputs
  integer, intent(in) :: v_dim1, v_dim2
  integer, intent(in) :: n1, n2
  real(8), intent(in) :: v(v_dim1, v_dim2)

  ! Inout arrays (will be reordered in-place)
  integer, intent(inout) :: inn_new1(n1), inn_new2(n2)

  ! Reorder both polygons
  call reorder_polygon(v, v_dim1, inn_new1, n1)
  call reorder_polygon(v, v_dim1, inn_new2, n2)

end subroutine ArrangeVertices


! -------------------------------------------------------
subroutine reorder_polygon(v, v_dim1, inn, n)
  implicit none
  integer, intent(in) :: v_dim1, n
  real(8), intent(in) :: v(v_dim1, *)
  integer, intent(inout) :: inn(n)

  integer :: i, j, idx_min
  real(8), allocatable :: x(:), y(:), ang(:)
  integer, allocatable :: order(:)

  real(8) :: cx, cy, dx, dy

  ! Temporary arrays
  allocate(x(n), y(n), ang(n), order(n))

  ! Extract coordinates in the current order
  do i = 1, n
     x(i) = v(1, inn(i))
     y(i) = v(2, inn(i))
  end do

  ! centroid
  cx = sum(x) / dble(n)
  cy = sum(y) / dble(n)

  ! angles relative to centroid (do not modify ang later)
  do i = 1, n
     ang(i) = atan2(y(i) - cy, x(i) - cx)
  end do

  ! initial order 1..n
  order = [(i, i = 1, n)]

  ! sort 'order' by ang(order(i)) descending => clockwise
  call sort_order_by_angle_desc(ang, order, n)

  ! apply sorted order to inn
  inn = inn(order)

  ! Now find index of lowest x (tie-breaker: lowest y) among the reordered list
  idx_min = 1
  do i = 2, n
     if ( v(1, inn(i)) < v(1, inn(idx_min)) ) then
        idx_min = i
     else if ( abs(v(1, inn(i)) - v(1, inn(idx_min))) < 1.0d-12 ) then
        ! tie in x -> choose lower y
        if ( v(2, inn(i)) < v(2, inn(idx_min)) ) idx_min = i
     end if
  end do

  ! rotate inn so that element at idx_min becomes first
  call rotate_array_to_k(inn, idx_min, n)

  deallocate(x, y, ang, order)
end subroutine reorder_polygon


! -------------------------------------------------------
subroutine sort_order_by_angle_desc(ang, order, n)
  implicit none
  real(8), intent(in) :: ang(n)
  integer, intent(inout) :: order(n)
  integer, intent(in) :: n
  integer :: i, j, tmp

  ! Simple O(n^2) sort swapping only 'order' based on ang(order)
  do i = 1, n - 1
     do j = i + 1, n
        if ( ang(order(i)) < ang(order(j)) ) then
           tmp = order(i)
           order(i) = order(j)
           order(j) = tmp
        end if
     end do
  end do
end subroutine sort_order_by_angle_desc


! -------------------------------------------------------
subroutine rotate_array_to_k(arr, k, n)
  implicit none
  integer, intent(inout) :: arr(n)
  integer, intent(in) :: k, n
  integer :: i
  integer, allocatable :: tmp(:)

  allocate(tmp(n))
  ! tmp(1) = arr(k), tmp(2) = arr(k+1), ... wrap around
  do i = 1, n
     tmp(i) = arr( mod(k-1 + (i-1), n) + 1 )
  end do
  arr = tmp
  deallocate(tmp)
end subroutine rotate_array_to_k



subroutine Squeeze_Tissue
  implicit none
  real*8 :: lowest_vy


  print*, 'Squeeze tissue ; percentage squeeze', percent_squeeze, '%'

  call Find_boundary_dynamic
  

  
  lowest_vy = minval((v(2,bottom_border(1:bottom_border_count))))

  v(2,:) = (v(2,:) - lowest_vy) * (1.0d0 - percent_squeeze/100.d0) + lowest_vy




end subroutine Squeeze_Tissue









end module Geometry
