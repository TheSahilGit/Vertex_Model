module T1_swap
  use array_info
  use Geometry
  use allocation

  integer, allocatable, dimension(:) :: Affected, Occurrences

  logical :: T1_pass


  contains


    subroutine find_T1
      implicit none


      cell_no = 0
      ver_no = 0
      d_val = 0


      count_T1 = 0
      do ic = 1, Nc !Lx*Ly

        if(num(ic).le.3)cycle


        call Gather_Cell_Vertices_PBC(inn(1:num(ic), ic), num(ic), vx, vy)

        do jc = 1, num(ic)

          next_idx = jc + 1
          prev_idx = jc - 1

          if(jc == num(ic))then
            next_idx = 1
          elseif(jc == 1)then
            prev_idx = num(ic)
          end if

          dx = vx(jc) - vx(next_idx)
          dy = vy(jc) - vy(next_idx)
          len_d_sq = dx*dx + dy*dy

          if(len_d_sq < min_d_T1*min_d_T1)then
            len_d = sqrt(len_d_sq)
            count_T1 = count_T1 + 1
            cell_no(count_T1) = ic
            ver_no(count_T1) = jc
            ver_no_next(count_T1) = next_idx
            d_val(count_T1) = len_d
          end if

        end do

      end do

      if(count_T1>0)then
        call random_number(rr)
        chosen_index = int(rr*count_T1 + 1)

         cellNoT1 = cell_no(chosen_index)
         verNoT1 = ver_no(chosen_index)
         verNoNextT1 = ver_no_next(chosen_index)
         len_d_T1 = d_val(chosen_index)
        
      !  cellNoT1 = cell_no(1)
      !  verNoT1 = ver_no(1)
      !  verNoNextT1 = ver_no_next(1)
      !  len_d_T1 = d_val(1)

      else
        cellNoT1 = 0
        verNoT1 = 0
        verNoNextT1 = 0
        len_d_T1 = 0.0d0
        T1_pass = .false.

      end if

!      write(*,*)'ver&cell', cellNoT1, verNoT1


      
   end subroutine find_T1

   subroutine find_T1_Affected

     implicit none

     integer :: ik, jcNext
     integer :: temp_array(Nc) !(Lx*Ly)
     

    ! write(*,*)inn(verNoT1,cellNoT1), inn(verNoNextT1, cellNoT1)
     

    ik = 0
     do ic = 1, Nc !Lx*Ly
       do jc = 1, num(ic)
         if((inn(jc,ic)==inn(verNoT1,cellNoT1)).or. & 
           inn(jc,ic)==inn(verNoNextT1,cellNoT1))then
        
       ! jcNext = jc + 1
       ! if(jc == num(ic))then
       !   jcNext = 1
       ! end if

       ! if((inn(jc,ic)==inn(verNoT1,cellNoT1)).and. &
       !    (inn(jcNext, ic) == inn(verNoNextT1,cellNoT1)) .or. &
       !    inn(jc,ic)==inn(verNoNextT1,cellNoT1) .and. &
       !    (inn(jcNext,ic)==inn(verNoT1,cellNoT1))) then

           ik = ik + 1
            
           temp_array(ik) = ic

         end if

       end do

     end do


     call find_unique_sorted(temp_array, ik, Affected, Occurrences)

!     write(*,*)size(Affected)
!     T1_pass = .true.
!     if(size(Affected).ne.4)then
!       T1_pass = .false.
!     end if


   end subroutine find_T1_Affected

   subroutine T1_core
     implicit none

     real*8 :: scaling_factor
     real*8 :: OldX(2), OldY(2), NewX(2), NewY(2)
     integer :: ik,il, cellNo, verNo_in1, verNo_in2, verNo_in1_indx, verNo_in2_indx
     integer :: inn_affected(inn_dim1), inn_temp(inn_dim1)

     inn_affected = 0
     


      !write(*,*)Affected
      !write(*,*)Occurrences


     scaling_factor = 2.0d0*min_d_T1/len_d_T1

!     write(*,*)cellNoT1, verNoT1, verNoNextT1, len_d_T1, scaling_factor


    OldX(1) = v(1, inn(verNoT1, cellNoT1))
    OldX(2) = v(1, inn(verNoNextT1, cellNoT1))
    OldY(1) = v(2, inn(verNoT1, cellNoT1))
    OldY(2) = v(2, inn(verNoNextT1, cellNoT1))

    ! PBC (log.txt): OldX(2)/OldY(2) must be in the same local unwrapped
    ! frame as OldX(1)/OldY(1) before the midpoint/bisector below, or a
    ! genuinely short edge that happens to straddle a periodic wrap would
    ! get a nonsensical (near-box-spanning) midpoint instead.
    call Unwrap_Relative(OldX(1), OldY(1), OldX(2), OldY(2))

    call PerpendicularBisector2(OldX,OldY,scaling_factor,NewX,NewY)

    ! PBC: NewX/NewY were computed in OldX(1)/OldY(1)'s local frame, which
    ! may itself sit outside the canonical box (or the flip may push a
    ! vertex just past a periodic edge) -- re-wrap before storing.
    call Wrap_Position(NewX(1), NewY(1))
    call Wrap_Position(NewX(2), NewY(2))

    v(1, inn(verNoT1, cellNoT1)) = NewX(2)
    v(1, inn(verNoNextT1, cellNoT1)) = NewX(1)
    v(2, inn(verNoT1, cellNoT1)) = NewY(2)
    v(2, inn(verNoNextT1, cellNoT1)) = NewY(1)
        
    verNo_in1 = inn(verNoT1, cellNoT1)
    verNo_in2 = inn(verNoNextT1, cellNoT1)
    
!    write(*,*)'verNo', verNo_in1, verNo_in2

    do ik = 1,size(Affected)

    
      if(Occurrences(ik).eq.2)then   ! These will lose.
        inn_temp = 0
        inn_affected = 0

        cellNo = Affected(ik)
        inn_affected(1:num(cellNo)) = inn(1:num(cellNo), cellNo)

!        write(*,*)'inn Affected', inn_affected(1:num(cellNo))

        do il = 1, num(cellNo)
          if((inn_affected(il)).eq.verNo_in1)then
            verNo_in1_indx = il
          end if
          if((inn_affected(il)).eq.verNo_in2)then
            verNo_in2_indx = il
          end if
        end do

        if(verNo_in1_indx.eq.1.and.verNo_in2_indx.eq.num(cellNo))then
          verNo_in1_indx = num(cellNo) + 1
        elseif(verNo_in2_indx.eq.1.and.verNo_in1_indx.eq.num(cellNo))then
          verNo_in2_indx = num(cellNo) + 1
        end if

       ! write(*,*)verNo_in1_indx, verNo_in2_indx

       if(verNo_in2_indx.lt.verNo_in1_indx)then
          
         inn_temp = 0

         do il = 1, verNo_in2_indx-1
            inn_temp(il) = inn_affected(il)
         end do


         do il = verNo_in2_indx, inn_dim1-1
           inn_temp(il) = inn_affected(il+1)
         end do

         do il = 1, num(cellNo)
           inn(il, cellNo) = inn_temp(il)
         end do

         num(cellNo) = num(cellNo) - 1

       elseif(verNo_in1_indx.lt.verNo_in2_indx)then
        

         inn_temp = 0

         do il = 1, verNo_in1_indx-1
            inn_temp(il) = inn_affected(il)
         end do

        
         do il = verNo_in1_indx, inn_dim1-1
           inn_temp(il) = inn_affected(il+1)
         end do

         do il = 1, num(cellNo)
           inn(il, cellNo) = inn_temp(il)
         end do

         num(cellNo) = num(cellNo) - 1

       end if
        

!       write(*,*)'inn updated', inn(1:num(cellNo), cellNo) 






      else if(Occurrences(ik).eq.1)then  ! These will gain. 
        
       inn_temp = 0
       inn_affected = 0

       verNo_in1_indx = 0
       verNo_in2_indx = 0

        cellNo = Affected(ik)
        inn_affected(1:num(cellNo)) = inn(1:num(cellNo), cellNo)

!       write(*,*)'inn_affected1',inn_affected

        do il = 1, num(cellNo)
          if((inn_affected(il)).eq.verNo_in1)then
             verNo_in1_indx  = il
           end if
           if((inn_affected(il)).eq.verNo_in2)then
             verNo_in2_indx  = il
           end if
        end do
       ! write(*,*)'index',verNo_in1_indx, verNo_in2_indx
        
       
       ! BUGFIX (log.txt): guard against overflowing inn's first dimension
       ! (inn_dim1) when a cell that is already at max vertex capacity would
       ! gain a vertex from this T1 flip. Previously there was no check here,
       ! so num(cellNo) could exceed inn_dim1 and the write-back loop below
       ! would write past inn's declared bound.
       if(verNo_in1_indx.ne.0.and.verNo_in2_indx.eq.0)then

         if(num(cellNo) .ge. inn_dim1)then
           write(*,*)'T1_core: cell', cellNo, 'already at max vertex capacity (inn_dim1) -- skipping vertex gain'
         else

         inn_temp = 0

         do il = 1, verNo_in1_indx
           inn_temp(il) = inn_affected(il)
         end do

         inn_temp(verNo_in1_indx+1) = verNo_in2

         do il = verNo_in1_indx + 2, inn_dim1
           inn_temp(il) = inn_affected(il-1)
         end do

         num(cellNo) = num(CellNo) + 1

         do il = 1, num(cellNo)
            inn(il, cellNo) = inn_temp(il)
         end do

         end if

       end if


       if(verNo_in2_indx.ne.0.and.verNo_in1_indx.eq.0)then

         if(num(cellNo) .ge. inn_dim1)then
           write(*,*)'T1_core: cell', cellNo, 'already at max vertex capacity (inn_dim1) -- skipping vertex gain'
         else

         inn_temp = 0
         do il = 1, verNo_in2_indx
           inn_temp(il) = inn_affected(il)
         end do

         inn_temp(verNo_in2_indx+1) = verNo_in1

         do il = verNo_in2_indx + 2, inn_dim1
           inn_temp(il) = inn_affected(il-1)
         end do

         num(cellNo) = num(CellNo) + 1

         do il = 1, num(cellNo)
            inn(il, cellNo) = inn_temp(il)
         end do

         end if

       end if
      
!       write(*,*)'inn_updated1',inn(:, cellNo)



     end if  ! Occurrence if
    
   end do   ! ik do
     



    end subroutine T1_core



    subroutine Do_T1


      implicit none
      integer :: il, im
      logical :: affected_computed

      ! Event-log locals (log.txt: T1/T2 spatial/cell-identity tracking, for
      ! Movie_Code.m to later overlay). Only used right before T1_core, once
      ! T1_pass has survived every rejection check above. ev_ids holds the
      ! numeric suffix of each affected cell's persistent cell_identity
      ! (CellIdNum, allocation.f90), not the string itself -- see
      ! Open_Event_Logs' comment for why.
      real*8 :: ev_x1, ev_y1, ev_x2, ev_y2, ev_mx, ev_my
      real*8 :: ev_ids(4)
      integer :: ev_ik

      T1_pass = .true.
      affected_computed = .false.

      call find_T1

      if(T1_pass)then
        if(if_Fixed_boundary)then
          call find_T1_Affected
          affected_computed = .true.


       if(sum(Total_T2_count(1:it)).gt.0)then
         call Find_boundary_dynamic
       else
         call Get_Boundary_info
       end if



          do im = 1, size(boundary)
            do il = 1, size(Affected)
              if(boundary(im).eq.Affected(il))then
                  T1_pass = .false.
                  write(*,*)'Boundary Ignored T1'
              end if
            end do
          end do
       end if
      end if

      
      if(T1_pass .and. (if_bottom_borders_fixed .or. if_top_borders_fixed)) then
        call Find_boundary_dynamic
      end if

      ! BUGFIX (log.txt): the second endpoint of the shrinking edge is
      ! verNoNextT1 (computed with wraparound in find_T1), not the raw index
      ! verNoT1+1. When the short edge is the polygon-closing edge
      ! (verNoT1 == num(cellNoT1)), inn(verNoT1+1, cellNoT1) read a stale/zero
      ! slot instead of the real second endpoint, letting a T1 flip through a
      ! fixed border undetected (and risking an out-of-bounds read on inn if
      ! the cell is at max vertex capacity).
      if(T1_pass.and.if_bottom_borders_fixed)then
      !!call Find_boundary_dynamic
        do im = 1, bottom_border_count
          if(bottom_border(im).eq.inn(verNoT1, cellNoT1) &
            .or.bottom_border(im).eq.inn(verNoNextT1, cellNoT1))then
            T1_pass = .false.
            write(*,*)'Bottom Border Ignored T1'
          end if
        end do
      end if

      if(T1_pass.and.if_top_borders_fixed)then
      !!call Find_boundary_dynamic
        do im = 1, top_border_count
          if(top_border(im).eq.inn(verNoT1, cellNoT1) &
            .or.top_border(im).eq.inn(verNoNextT1, cellNoT1))then
            T1_pass = .false.
            write(*,*)'Top Border Ignored T1'
          end if
        end do
      end if

     
!      write(*,*)T1_pass,count_T1
      
      if(T1_pass)then
          ! OPTIMIZATION (log.txt): find_T1_Affected (an O(Nc*max_verts)
          ! full-mesh scan) was called unconditionally here even when the
          ! if_Fixed_boundary branch above already called it with the same
          ! inputs (cellNoT1/verNoT1/verNoNextT1 unchanged since), doubling
          ! the scan for every accepted T1. Skip it if already computed.
          if (.not. affected_computed) call find_T1_Affected

          ! BUGFIX (log.txt): reject a T1 flip that would reduce any affected
          ! "losing" cell (Occurrences==2, i.e. it shares the collapsing
          ! edge) below the minimum viable polygon vertex count (3).
          ! T1_core's Occurrences==2 branch had no lower-bound guard -- a
          ! cell already at 3 vertices could be flipped down to 2, and later
          ! down to 1, at which point num(ic)==1 makes prev_idx wrap to 0 in
          ! Force_Calculation (an out-of-bounds array access). Confirmed live
          ! under an aggressive min_d_T1 stress test. A cell this small
          ! should be removed by T2, not further reduced by T1 -- reject the
          ! whole flip (no partial mutation) rather than trying to patch just
          ! this cell's update.
          do il = 1, size(Affected)
            if (Occurrences(il) == 2 .and. num(Affected(il)) <= 3) then
              T1_pass = .false.
              write(*,*)'Degenerate cell Ignored T1'
            end if
          end do

          if (T1_pass) then
            ! Event log (log.txt): capture the flipping edge's midpoint and
            ! the (up to 4) affected cells' PERSISTENT identities *before*
            ! T1_core overwrites the edge's vertex positions in place --
            ! afterward, v(:,inn(verNoT1/verNoNextT1,cellNoT1)) holds the
            ! POST-flip positions, not the event location. PBC-unwrapped
            ! (Unwrap_Relative/Wrap_Position), same convention as every
            ! other geometric routine in this codebase.
            ev_x1 = v(1, inn(verNoT1,     cellNoT1)); ev_y1 = v(2, inn(verNoT1,     cellNoT1))
            ev_x2 = v(1, inn(verNoNextT1, cellNoT1)); ev_y2 = v(2, inn(verNoNextT1, cellNoT1))
            call Unwrap_Relative(ev_x1, ev_y1, ev_x2, ev_y2)
            ev_mx = 0.5d0 * (ev_x1 + ev_x2)
            ev_my = 0.5d0 * (ev_y1 + ev_y2)
            call Wrap_Position(ev_mx, ev_my)

            ev_ids = 0.0d0
            do ev_ik = 1, min(4, size(Affected))
              ev_ids(ev_ik) = CellIdNum(cell_identity(Affected(ev_ik)))
            end do
            write(iunit_T1events) dble(it), ev_mx, ev_my, &
              ev_ids(1), ev_ids(2), ev_ids(3), ev_ids(4)

            call T1_core
           ! print*, "T1 Happened at it = ", it
            Total_T1_count(it) = Total_T1_count(it) + 1
          end if
        end if



    end subroutine Do_T1






end module T1_swap
