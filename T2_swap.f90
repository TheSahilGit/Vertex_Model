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


      cell_no = 0
      ver_no = 0
      area_val = 0

      count_T2 = 0

      do ic = 1, Lx*Ly
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
      integer :: temp_array(Lx*Ly)

      integer :: ik

      !write(*,*)cellNoT2

      ik = 0
      do ic = 1,Lx*Ly
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
      integer :: whole_inn_array_temp(inn_dim1, inn_dim2)
      integer :: whole_num_array_temp(num_dim)

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

    

    whole_inn_array_temp = 0
    whole_inn_array_temp(:,1:cellNoT2-1) = inn(:,1:cellNoT2-1)
    whole_inn_array_temp(:, cellNoT2:inn_dim2-1) = inn(:,cellNoT2+1:inn_dim2)

    whole_num_array_temp = 0
    whole_num_array_temp(1:cellNoT2-1)  = num(1:cellNoT2-1)
    whole_num_array_temp(cellNoT2 : num_dim-1) = num(cellNoT2+1:num_dim)
   
    inn = whole_inn_array_temp
    num = whole_num_array_temp


    end subroutine T2_core

   subroutine Do_T2

     implicit none

     integer :: im, il

     T2_pass = .true.



    call find_T2

    if(T2_pass.and.if_Fixed_boundary)then

      call find_T2_Affected
      call Get_Boundary_info


      do im = 1, size(boundary)
        do il = 1, size(Affected_T2)
          if(boundary(im).eq.Affected_T2(il))then
            T2_pass = .false.
            write(*,*)"Boundary Ignored T2"
          end if
        end do
      end do

    end if
    
    if(T2_pass)then
       call find_T2_Affected
       call T2_core
       Total_T2_count = Total_T2_count + 1
     end if



   
   end subroutine Do_T2




end module T2_swap
