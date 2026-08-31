module T2_swap
  use array_info
  use Geometry
  use allocation


  integer, allocatable, dimension(:) :: Affected_T2, Occurrences_T2
  integer :: index_array(128) ! 128 has no reason. 

  logical :: T2_pass

  contains

    subroutine find_T2
      implicit none
      integer :: im, il


      cell_no = 0
      ver_no = 0
      area_val = 0

      count_T2 = 0

      do ic = 1, Nc !Lx*Ly
        if(num(ic).eq.3)then
          vx = v(1,inn(1:num(ic),ic))
          vy = v(2,inn(1:num(ic),ic))
          call  CalculateArea(vx,vy,num(ic),area)
          area = abs(area)

          if (area.le.min_area_T2)then
            count_T2 = count_T2 + 1
            area_val(count_T2) = area
            cell_no(count_T2) = ic
          end if

        end if

      end do


      if(count_T2>0)then
        call random_number(rr)
        chosen_index = int(rr*count_T2+1)
        cellNoT2 = cell_no(chosen_index)
      else
        T2_pass = .false.
      end if


    end subroutine find_T2
    

    subroutine find_T2_Affected
      implicit none
      integer :: temp_array(Nc) !(Lx*Ly)

      integer :: ik

      !write(*,*)cellNoT2

      ik = 0
      do ic = 1, Nc !Lx*Ly
        do jc = 1, num(ic)

         ! write(*,*)inn(1,cellnoT2)
          if(inn(jc,ic) == inn(1,cellnoT2).or. &
            inn(jc,ic) == inn(2,cellnoT2) .or. &
            inn(jc,ic) == inn(3,cellnoT2)) then
            ik = ik + 1
            temp_array(ik) = ic
            index_array(ik) = jc
          end if

        end do
      end do


      call find_unique_sorted(temp_array, ik, Affected_T2, Occurrences_T2)

!      write(*,*)'affected_cell',cellNoT2
!      write(*,*)'all_affected_cell',Affected_T2
!      write(*,*)'occurrence', Occurrences_T2
!      write(*,*)'index_array', index_array(1:ik)
      


    end subroutine find_T2_Affected

    subroutine T2_core

      implicit none

      integer :: ik, il, im

      integer ::  cellNo_aff

      real*8 :: VcmX, VcmY
      integer :: inn_affected(inn_dim1),inn_temp(inn_dim1)

      integer :: start_index, stop_index

      vx = v(1,inn(1:num(cellNoT2), cellNoT2))
      vy = v(2,inn(1:num(cellNoT2), cellNoT2))

      VcmX = sum(vx)/size(vx)
      VcmY = sum(vy)/size(vy)

      ! Replacing 1st vertex of cellNoT2 with the COM,
      ! and will update this and delete rest. 

      v(1, inn(1,cellNoT2)) = VcmX
      v(2, inn(1,cellNoT2)) = VcmY

      verNoT2 =  inn(1,cellNoT2) ! This goes in the rest of the cells

      ! Occurrences = 3 --> The cell itself. 
      ! Occurrences = 2 --> Shares an edge.
      ! Occurrences = 1 --> Shares a vertex.

      start_index = 1

      do ik = 1,size(Affected_T2)


        inn_affected = 0
        inn_temp = 0
        
        cellNo_aff = Affected_T2(ik)
        
        stop_index = start_index + Occurrences_T2(ik) - 1

        if(Occurrences_T2(ik).eq.3) goto 31

        inn_affected(1:num(cellNo_aff)) = & 
          inn(1:num(cellNo_aff), cellNo_aff)


        inn_temp(1:index_array(start_index)-1) =  & 
          inn_affected(1:index_array(start_index)-1)

        inn_temp(index_array(start_index)) = verNoT2  

        if(index_array(stop_index).eq.num(cellNo_aff))then
          inn_temp(index_array(start_index)+1 : num(cellNo_aff) - & 
             Occurrences_T2(ik)+2) = inn_affected(index_array(start_index)+1 : &
             num(cellNo_aff))

        else
          inn_temp(index_array(start_index)+1 : num(cellNo_aff) - & 
            Occurrences_T2(ik)+1) =  inn_affected(index_array(stop_index)+1 :  &
            num(cellNo_aff))
           
         end if

        num(cellNo_aff) = num(cellNo_aff) - Occurrences_T2(ik)+ 1 

        inn(1:num(cellNo_aff), cellNo_aff) = inn_temp(1:num(cellNo_aff))

        inn(num(cellNo_aff)+1:inn_dim1, cellNo_aff) = 0
        

31        start_index = stop_index+1
    

    end do

    

    ! OPTIMIZATION (log.txt): previously this rebuilt a whole_inn_array_temp
    ! (inn_dim1 x inn_dim2) / whole_num_array_temp/Rho_temp/ROCK_temp/
    ! Myosin_temp/cell_identity_temp (num_dim, incl. a character(500) array)
    ! from scratch every T2 event (up to 500x/timestep), even though only
    ! cells 1:Nc are ever live -- columns Nc+1:inn_dim2/num_dim are always
    ! zero/'cell_0' and stay that way. Shifting in-place over 1:Nc only
    ! (Fortran array assignment evaluates the RHS before assigning, so this
    ! overlapping shift is well-defined) gives an identical result at a
    ! fraction of the cost, with no large temporaries.
    inn(:, cellNoT2:Nc-1) = inn(:, cellNoT2+1:Nc)
    inn(:, Nc) = 0

    num(cellNoT2:Nc-1) = num(cellNoT2+1:Nc)
    num(Nc) = 0


    if(if_RhoROCK)then

      Rho(cellNoT2:Nc-1) = Rho(cellNoT2+1:Nc)
      Rho(Nc) = 0.0d0

      ROCK(cellNoT2:Nc-1) = ROCK(cellNoT2+1:Nc)
      ROCK(Nc) = 0.0d0

      Myosin(cellNoT2:Nc-1) = Myosin(cellNoT2+1:Nc)
      Myosin(Nc) = 0.0d0

    end if


    cell_identity(cellNoT2:Nc-1) = cell_identity(cellNoT2+1:Nc)
    cell_identity(Nc) = 'cell_0'


    end subroutine T2_core

   subroutine Do_T2
     implicit none

     integer :: im, il
     logical :: affected_computed

     T2_pass = .true.
     affected_computed = .false.



    call find_T2

     if(T2_pass.and.if_Fixed_boundary)then

       call find_T2_Affected
       affected_computed = .true.

      if(sum(Total_T2_count(1:it)).gt.0)then
        call Find_boundary_dynamic
      else
        call Get_Boundary_info
      end if

 
 
       do im = 1, size(boundary)
         do il = 1, size(Affected_T2)
           if(boundary(im).eq.Affected_T2(il))then
             T2_pass = .false.
             write(*,*)"Boundary Ignored T2"
           end if
         end do
       end do

     end if

     ! BUGFIX (log.txt): Do_T1 explicitly protects bottom_border/top_border
     ! pinned vertices (if_bottom_borders_fixed/if_top_borders_fixed) in
     ! addition to the generic if_Fixed_boundary check above; Do_T2 had no
     ! equivalent, so a T2 extrusion could remove a cell touching a "fixed"
     ! border and silently break that boundary condition.
     if(T2_pass .and. (if_bottom_borders_fixed .or. if_top_borders_fixed)) then
       call Find_boundary_dynamic
     end if

     if(T2_pass.and.if_bottom_borders_fixed)then
       do im = 1, bottom_border_count
         if(bottom_border(im).eq.inn(1,cellNoT2) .or. &
            bottom_border(im).eq.inn(2,cellNoT2) .or. &
            bottom_border(im).eq.inn(3,cellNoT2))then
           T2_pass = .false.
           write(*,*)'Bottom Border Ignored T2'
         end if
       end do
     end if

     if(T2_pass.and.if_top_borders_fixed)then
       do im = 1, top_border_count
         if(top_border(im).eq.inn(1,cellNoT2) .or. &
            top_border(im).eq.inn(2,cellNoT2) .or. &
            top_border(im).eq.inn(3,cellNoT2))then
           T2_pass = .false.
           write(*,*)'Top Border Ignored T2'
         end if
       end do
     end if


    if(T2_pass)then
       ! OPTIMIZATION (log.txt): find_T2_Affected (an O(Nc*max_verts)
       ! full-mesh scan) was called unconditionally here even when the
       ! if_Fixed_boundary branch above already called it with the same
       ! inputs (cellNoT2 unchanged since), doubling the scan for every
       ! accepted T2. Skip it if already computed.
       if (.not. affected_computed) call find_T2_Affected

       ! BUGFIX (log.txt, re-review pass): mirrors the T1_core degenerate-
       ! cell guard added earlier for T1_swap.f90. T2_core's Occurrences_T2
       ! == 2 case (a neighbor sharing an edge with the collapsing triangle)
       ! decrements that neighbor's vertex count by 1 with no lower-bound
       ! check. A neighbor already at 3 vertices would be pushed to 2, and a
       ! second such T2 event touching the same cell later could push it to
       ! 1 -- at num(ic)==1, Force_Calculation's wraparound leaves
       ! prev_idx==0, an out-of-bounds array access. Reject the whole T2
       ! event instead (no partial mutation), consistent with the T1 fix.
       do il = 1, size(Affected_T2)
         if (Occurrences_T2(il) == 2 .and. num(Affected_T2(il)) <= 3) then
           T2_pass = .false.
           write(*,*)'Degenerate cell Ignored T2'
         end if
       end do

       if (T2_pass) then
         call T2_core
         Nc = Nc - 1
         Total_T2_count(it) = Total_T2_count(it) + 1
         print*, 'T2_happened, Nc = ', Nc
       end if
     end if



   
   end subroutine Do_T2




end module T2_swap
