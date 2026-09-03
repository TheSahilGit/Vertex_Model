module Proliferation
  use allocation
  use Geometry


  contains 


    subroutine Do_Proliferation

      implicit none

      !real*8 :: area_0
      real*8 :: prin_axis_x1, prin_axis_y1
      real*8 :: prin_axis_x2, prin_axis_y2
      real*8 :: x_intersection(2), y_intersection(2)
      integer :: n_intersections
      integer :: idx_pair(2,2)
      integer :: ic, icand
      logical :: division_success


      !area_0 = 1.25d0*Ao

      call  Find_Proliferation(area_0, v, inn, num, Nc, chosen_cell, chosen_cell_count)


     ! print*, 'chosen cell', chosen_cell(1:chosen_cell_count)

      ! BUGFIX (log.txt): previously only chosen_cell(1) was ever attempted.
      ! If that specific cell's division always fails (e.g. it sits on the
      ! tissue boundary and one of its split edges has no neighboring cell,
      ! see the othercells_count guard in Proliferation_Core), its area
      ! never shrinks, so it stays chosen_cell(1) on every subsequent
      ! timestep -- permanently starving every other over-threshold cell of
      ! ever being tried (confirmed live: cell division stopped completely,
      ! 0/2000 successes, once if_Fixed_boundary made the very first
      ! candidate perpetually un-dividable). Fall through to the next
      ! candidate whenever one fails, so a single un-dividable cell can no
      ! longer block the rest.
      if(chosen_cell_count.gt.0)then

        do icand = 1, chosen_cell_count

        ic = chosen_cell(icand)

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

        division_success = .false.

        ! BUGFIX (log.txt): n_intersections was computed but never checked.
        ! If the bisector doesn't cleanly cross exactly 2 edges (degenerate
        ! bisector, or it passes through a vertex), idx_pair is left partly
        ! uninitialized and Proliferation_Core/SplitPolygon can end up zeroing
        ! out the dividing cell's vertex count. Skip the division instead.
        if(n_intersections .ne. 2)then
          write(*,*)'Proliferation skipped: bisector did not cleanly intersect 2 edges for cell', &
            ic, ' (n_intersections =', n_intersections, ')'
        else

         call Proliferation_Core(ic, Nc, v, inn, num, &
              num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2, &
              x_intersection, y_intersection, idx_pair, division_success)

         !print*, num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2

        end if

        if (division_success) exit

        end do

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
               x_intersection, y_intersection, idx_pair, success)

      implicit none

      integer, intent(in) :: ic
      real*8, intent(in) :: x_intersection(2), y_intersection(2)
      integer, intent(in) :: num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2
      integer, intent(in) :: idx_pair(2,2)
      integer, intent(inout) :: Nc
      real*8, intent(inout) :: v(v_dim1, v_dim2)
      integer, intent(inout) :: inn(inn_dim1, inn_dim2)
      integer, intent(inout) :: num(num_dim)
      logical, intent(out) :: success

      integer :: maxinn, ii, jc
      integer :: innaff_pair1(2), innaff_pair2(2)


      integer :: othercells(num_dim), othercells_count
      integer :: inn_new1(num(ic)+2), inn_new2(num(ic)+2)
      integer :: n1, n2, n1_new, n2_new

      integer :: inn_neighbor1(inn_dim1), inn_neighbor2(inn_dim1)
      ! BUGFIX (log.txt): UpdateNeighborPolygons (Geometry.f90) declares its
      ! output dummy args as size(inn_old1)+2 / size(inn_old2)+2 -- and
      ! inn_neighbor1/inn_neighbor2 are passed to it WHOLE (unsliced), so
      ! size(inn_old1) inside that subroutine is inn_dim1, not
      ! num(othercells(1)). It therefore always writes up to inn_dim1+2
      ! elements into inn_neighbor1_new/inn_neighbor2_new, every single
      ! division, regardless of how many vertices the neighbor cell
      ! actually has. These were only sized inn_dim1 (no headroom at all),
      ! so every division overflowed them by exactly 2 integers -- caught
      ! live via AddressSanitizer (heap-buffer-overflow WRITE at
      ! Geometry.f90:642, `inn_new1 = 0`) once a large-mesh run made the
      ! corrupted heap bytes matter enough to crash on a later, unrelated
      ! free() ("free(): invalid pointer"). The n1_new/n2_new > inn_dim1
      ! guard further down only rejects the *result* after this overflow
      ! had already happened -- it could never prevent the underlying
      ! out-of-bounds write. Sized with the same +2 headroom the callee
      ! actually needs.
      integer :: inn_neighbor1_new(inn_dim1+2), inn_neighbor2_new(inn_dim1+2)
      




      !print*, "***********"
      !print*, Nc, num_dim, inn_dim1, inn_dim2, v_dim1, v_dim2
      !print*,  x_intersection, y_intersection

      success = .false.

      maxinn =  maxval(inn)

      !print*, 'maxinn', maxinn

      ! BUGFIX (log.txt): no capacity check previously existed before growing
      ! the vertex/cell arrays here -- once pre-allocated headroom (v_dim2 for
      ! vertices, inn_dim2/num_dim for cells) was exhausted, this silently
      ! wrote out of bounds. Bail out cleanly instead.
      if (maxinn+2 .gt. v_dim2) then
        write(*,*)'Proliferation_Core: skipped -- dividing cell', ic, &
          'would exceed v_dim2 vertex capacity (maxinn+2 =', maxinn+2, ', v_dim2 =', v_dim2, ')'
        return
      end if
      if (Nc+1 .gt. inn_dim2 .or. Nc+1 .gt. num_dim) then
        write(*,*)'Proliferation_Core: skipped -- dividing cell', ic, &
          'would exceed inn_dim2/num_dim cell capacity (Nc+1 =', Nc+1, ')'
        return
      end if

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

        ! BUGFIX (log.txt): the code below unconditionally assumes exactly 2
        ! neighboring cells were found (othercells(1), othercells(2)). If one
        ! of the dividing cell's two split edges lies on the tissue's outer
        ! boundary (a realistic case for a boundary cell), othercells_count
        ! can be 0 or 1, leaving othercells(2) (or (1)) == 0 -- an invalid
        ! index-0 access into inn/num just below. Bail out cleanly instead.
        if (othercells_count .ne. 2) then
          write(*,*)'Proliferation_Core: skipped -- cell', ic, &
            'has a split edge on the tissue boundary (othercells_count =', othercells_count, ')'
          return
        end if


      call SplitPolygon(inn(1:num(ic), ic), num(ic), idx_pair(1,:), idx_pair(2,:), &
        maxinn+1, maxinn+2, inn_new1, n1, inn_new2, n2)

      ! BUGFIX (log.txt, re-review pass): SplitPolygon has its own internal
      ! failure paths (chord vertices not found in inn_ic, or an "overflow
      ! constructing inn_new1/inn_new2" guard) that leave n1/n2 invalid (0,
      ! or a partial count above num(ic)+2) and return without completing
      ! the split. Without this check, execution fell through to
      ! ArrangeVertices (which would then read inn_new1(1:n1) out of its
      ! declared bound in the overflow case) and eventually
      ! inn(1:n1,ic)=inn_new1(1:n1)/num(ic)=n1, corrupting the dividing
      ! cell's topology while this subroutine still reported success=.true.
      ! at the end. Bail out cleanly instead, exactly like the other guards.
      if (n1 < 3 .or. n1 > num(ic)+2 .or. n2 < 3 .or. n2 > num(ic)+2) then
        write(*,*)'Proliferation_Core: skipped -- SplitPolygon failed for cell', ic, &
          '(n1 =', n1, ', n2 =', n2, ')'
        return
      end if

!      print*, 'newinn1', n1, '-', inn_new1
!      print*, 'newinn2', n2, '-', inn_new2


      call ArrangeVertices(v, v_dim1, v_dim2, inn_new1, n1, inn_new2, n2)


       !print*, othercells(1:othercells_count)

       inn_neighbor1(1:num(othercells(1))) = inn(1:num(othercells(1)), othercells(1))
       inn_neighbor2(1:num(othercells(2))) = inn(1:num(othercells(2)), othercells(2))


!       print*, 'neighbor1',  inn_neighbor1(1:num(othercells(1)))
!       print*, 'neighbor2', inn_neighbor2(1:num(othercells(2)))
!
!       print*,'pair1',  idx_pair(1,:)
!       print*,'pair2',  idx_pair(2,:)
!       print*, maxinn+1, maxinn+2
       

      ! inn_neighbor1/2 are passed whole (not sliced to num(othercells(.))),
      ! so UpdateNeighborPolygons' output dummies are sized inn_dim1+2 --
      ! the actual arguments below must match that exactly (see the
      ! inn_neighbor1_new/inn_neighbor2_new declaration comment above).
      call UpdateNeighborPolygons(inn_neighbor1, num(othercells(1)), inn_neighbor2, num(othercells(2)), &
                                   idx_pair(1,:), idx_pair(2,:), maxinn+1, maxinn+2, &
                                   inn_neighbor1_new, n1_new, &
                                   inn_neighbor2_new, n2_new)

      ! BUGFIX (log.txt, re-review pass): no check previously existed that
      ! n1/n2 (the dividing cell's two daughters) or n1_new/n2_new (the two
      ! neighboring cells, each grown by 1 vertex) stay within inn_dim1 --
      ! the hard cap on vertices per cell (inn(inn_dim1, inn_dim2)). A cell
      ! that reached inn_dim1 vertices via sustained T1 activity and later
      ! also exceeds the division-area threshold would silently overflow
      ! inn's first dimension on the writes below. Bail out cleanly instead,
      ! mirroring the analogous guard already added for T1_core.
      if (n1 > inn_dim1 .or. n2 > inn_dim1 .or. n1_new > inn_dim1 .or. n2_new > inn_dim1) then
        write(*,*)'Proliferation_Core: skipped -- cell', ic, &
          'division would exceed inn_dim1 vertex-per-cell capacity', &
          '(n1,n2,n1_new,n2_new =', n1, n2, n1_new, n2_new, ')'
        return
      end if

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


      ! Assigning motility to the two brand-new vertices created by this
      ! division (maxinn+1, on the pair1(1)-pair1(2) edge; maxinn+2, on the
      ! pair2(1)-pair2(2) edge -- see SplitPolygon). Every other vertex in
      ! both daughter cells is an already-existing mother-cell vertex, so it
      ! keeps whatever mot value it already had -- only these two are new
      ! and need one assigned. BUGFIX (log.txt): previously both were copied
      ! from mot(inn(1,ic)), a single arbitrary vertex of daughter cell ic
      ! unrelated to either split edge -- correct only when mot is spatially
      ! uniform (plain if_motility) and silently wrong under
      ! if_motility_gradient/if_motility_hotspot. Interpolate each new
      ! vertex from the actual two mother-cell vertices whose edge it splits
      ! instead, so both daughter cells' motility remains consistent with
      ! the mother's local field regardless of how it varies in space.
      mot(maxinn + 1) = 0.5d0 * (mot(idx_pair(1,1)) + mot(idx_pair(1,2)))
      mot(maxinn + 2) = 0.5d0 * (mot(idx_pair(2,1)) + mot(idx_pair(2,2)))


!      print*, 'after inn1', inn(1:num(ic), ic)
!      print*, 'after inn2', inn(1:num(Nc+1), Nc+1)
!      
!      print*, 'after v1', v(1,inn(1:num(ic), ic))
!      print*, 'after v2', v(1,inn(1:num(Nc+1), Nc+1))

      Nc = Nc + 1
      success = .true.

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
