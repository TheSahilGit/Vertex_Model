module Proliferation
  use allocation
  use Geometry


  contains 


    subroutine Do_Proliferation

      implicit none

      real*8 :: area_0
      real*8 :: prin_axis_x1, prin_axis_y1
      real*8 :: prin_axis_x2, prin_axis_y2
      real*8 :: x_intersection(2), y_intersection(2)
      integer :: n_intersections
      integer :: idx_pair(2,2)
      integer :: ic
      

      area_0 = 1.3d0*Ao

      call  Find_Proliferation(area_0, v, inn, num, Nc, chosen_cell, chosen_cell_count)


     ! print*, 'chosen cell', chosen_cell(1:chosen_cell_count)

      if(chosen_cell_count.gt.0)then 


        ic = chosen_cell(1)
  
        nn = num(ic)
        
    
        vx = v(1,inn(1:nn,ic))
        vy = v(2,inn(1:nn,ic))
  
  
        call Find_Principle_Axis(vx, vy, nn, &
          prin_axis_x1, prin_axis_y1, prin_axis_x2, prin_axis_y2)
  
       !print*, 'principle_axis', prin_axis_x1, prin_axis_y1, prin_axis_x2, prin_axis_y2
  
        call Find_Bisector_Intersections(vx, vy, nn, & 
             prin_axis_x1, prin_axis_y1, prin_axis_x2, prin_axis_y2, &
             x_intersection, y_intersection, n_intersections, inn(1:nn, ic), &
             ic, idx_pair)
  
        ! print*,'pb', n_intersections, x_intersection, y_intersection
  
  


         call Proliferation_Core(ic, Nc, v, inn, num, & 
              num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2, &
              x_intersection, y_intersection, idx_pair)
  
         !print*, num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2

       end if
  


    end subroutine Do_Proliferation
    






    subroutine Find_Proliferation(area0, v, inn, num, Nc, chosen_cell, chosen_cell_count)

      implicit none

      real*8, intent(in) :: area0
      real*8, intent(in) :: v(v_dim1, v_dim2)
      integer, intent(in) :: inn(inn_dim1, inn_dim2)
      integer, intent(in) :: num(num_dim)
      integer, intent(in) :: Nc

      integer, intent(out) :: chosen_cell(v_dim2)
      integer, intent(out) :: chosen_cell_count


      integer :: ic



      chosen_cell = 0
      chosen_cell_count = 0 

      do ic = 1, Nc !Lx*Ly
  
        nn = num(ic)
  
        vx = v(1,inn(1:nn,ic))
        vy = v(2,inn(1:nn,ic))
  
        call CalculateArea(vx,vy,nn,area)
        area = abs(area)
        if(area.gt.area0) then
          chosen_cell_count = chosen_cell_count + 1
          chosen_cell(chosen_cell_count)  = ic
        end if


      end do


    end subroutine Find_Proliferation


    subroutine Proliferation_Core(ic, Nc, v, inn, num, & 
               num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2, &
               x_intersection, y_intersection, idx_pair)

      implicit none

      integer, intent(in) :: ic
      real*8, intent(in) :: x_intersection(2), y_intersection(2)
      integer, intent(in) :: num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2
      integer, intent(in) :: idx_pair(2,2)
      integer, intent(inout) :: Nc
      real*8, intent(inout) :: v(v_dim1, v_dim2)
      integer, intent(inout) :: inn(inn_dim1, inn_dim2)
      integer, intent(inout) :: num(num_dim)

      integer :: maxinn, ii, jc
      integer :: innaff_pair1(2), innaff_pair2(2)


      integer :: othercells(num_dim), othercells_count
      integer :: inn_new1(num(ic)+2), inn_new2(num(ic)+2)
      integer :: n1, n2, n1_new, n2_new

      integer :: inn_neighbor1(inn_dim1), inn_neighbor2(inn_dim1)
      integer :: inn_neighbor1_new(inn_dim1), inn_neighbor2_new(inn_dim1)
      




      !print*, "***********"
      !print*, Nc, num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2
      !print*,  x_intersection, y_intersection

      maxinn =  maxval(inn)

      !print*, 'maxinn', maxinn

      v(1, maxinn+1) = x_intersection(1)
      v(2, maxinn+1) = y_intersection(1)

      v(1, maxinn+2) = x_intersection(2)
      v(2, maxinn+2) = y_intersection(2)

      !print*, 'inn(ic)', num(ic), inn(1:num(ic),ic) 
      !print*, 'inn(Nc+1)', inn(:,Nc+1) 

      !print*, 'idx_pair', idx_pair(1,:)
      !print*, 'idx_pair', idx_pair(2,:)







        call Find_Others_affected(ic, Nc, inn, num, idx_pair, othercells, othercells_count)

!        print*, 'othercells= ',othercells(1:othercells_count)
!        print*, 'ic =', ic


      call SplitPolygon(inn(1:num(ic), ic), num(ic), idx_pair(1,:), idx_pair(2,:), &
        maxinn+1, maxinn+2, inn_new1, n1, inn_new2, n2)


!      print*, 'newinn1', n1, '-', inn_new1
!      print*, 'newinn2', n2, '-', inn_new2


      call ArrangeVertices(v, v_dim1, v_dim2, inn_new1, n1, inn_new2, n2)
    


       inn_neighbor1(1:num(othercells(1))) = inn(1:num(othercells(1)), othercells(1))
       inn_neighbor2(1:num(othercells(2))) = inn(1:num(othercells(2)), othercells(2))


!       print*, 'neighbor1',  inn_neighbor1(1:num(othercells(1)))
!       print*, 'neighbor2', inn_neighbor2(1:num(othercells(2)))
!
!       print*,'pair1',  idx_pair(1,:)
!       print*,'pair2',  idx_pair(2,:)
!       print*, maxinn+1, maxinn+2
       

      call UpdateNeighborPolygons(inn_neighbor1, num(othercells(1)), inn_neighbor2, num(othercells(2)), &
                                   idx_pair(1,:), idx_pair(2,:), maxinn+1, maxinn+2, &
                                   inn_neighbor1_new(1:num(othercells(1))), n1_new, & 
                                   inn_neighbor2_new(1:num(othercells(2))), n2_new)



!      print*, 'After'
!      print*, inn_neighbor1_new(1:n1_new)
!      print*, inn_neighbor2_new(1:n2_new)
      

      call ArrangeVertices(v, v_dim1, v_dim2, & 
        inn_neighbor1_new(1:n1_new), n1_new, inn_neighbor2_new(1:n2_new), n2_new)



!      print*, 'After2'
!      print*, inn_neighbor1_new(1:n1_new)
!      print*, inn_neighbor2_new(1:n2_new)

      inn(1:n1, ic) = inn_new1(1:n1)
      inn(1:n2, Nc+1) = inn_new2(1:n2)
      inn(1:n1_new, othercells(1)) = inn_neighbor1_new(1:n1_new)
      inn(1:n2_new, othercells(2)) = inn_neighbor2_new(1:n2_new)
      num(ic) = n1
      num(othercells(1)) = n1_new
      num(othercells(2)) = n2_new
      num(Nc+1) = n2


!      print*, 'after inn1', inn(1:num(ic), ic)
!      print*, 'after inn2', inn(1:num(Nc+1), Nc+1)
!      
!      print*, 'after v1', v(1,inn(1:num(ic), ic))
!      print*, 'after v2', v(1,inn(1:num(Nc+1), Nc+1))

      Nc = Nc + 1

      print*, 'Proliferation: ic,  Nc = ', ic, Nc

    


      

      

    end subroutine Proliferation_Core



subroutine Find_Others_affected(ic, Nc, inn, num, idx_pair, othercells, othercells_count)
  implicit none
  integer, intent(in)  :: ic, Nc
  integer, intent(in)  :: idx_pair(2,2)
  integer, intent(in)  :: inn(inn_dim1, inn_dim2)
  integer, intent(in)  :: num(num_dim)
  integer, intent(out) :: othercells(num_dim)
  integer, intent(out) :: othercells_count

  integer :: i, j, n_verts
  integer :: v1, v2, v3, v4
  integer :: a, b
  logical :: found

  v1 = idx_pair(1,1)
  v2 = idx_pair(1,2)
  v3 = idx_pair(2,1)
  v4 = idx_pair(2,2)

  othercells_count = 0
  othercells = 0

  ! Loop over all cells
  do i = 1, Nc
     if (i == ic) cycle  ! skip the current cell
     n_verts = num(i)
     found = .false.

     ! Loop over edges of the cell
     do j = 1, n_verts
        a = inn(j, i)
        if (j == n_verts) then
           b = inn(1, i)   ! close the polygon
        else
           b = inn(j+1, i)
        end if

        ! Check if edge matches (v1,v2) or (v2,v1)
        if ( (a == v1 .and. b == v2) .or. (a == v2 .and. b == v1) ) then
           found = .true.
           exit
        end if

        ! Check if edge matches (v3,v4) or (v4,v3)
        if ( (a == v3 .and. b == v4) .or. (a == v4 .and. b == v3) ) then
           found = .true.
           exit
        end if
     end do

     if (found) then
        othercells_count = othercells_count + 1
        othercells(othercells_count) = i
     end if
  end do

end subroutine Find_Others_affected




end module Proliferation
