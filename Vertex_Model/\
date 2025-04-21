module T1_swap_new
  use array_info
  use Geometry
  use allocation

  integer :: Affected(4)

  contains


    subroutine find_T1_2
      implicit none

      integer :: next_idx, prev_idx, chosen_index
      real*8 :: len_d, dx, dy, rr
      integer*4, dimension(:), allocatable :: cell_no, ver_no, ver_no_next
      real*8, dimension(:), allocatable :: d_val

      allocate(d_val(Lx*Ly*6)) ! If all the cells have 6 sides on average,this is the max number of vertices possible with repetation.

      allocate(cell_no(maxval(inn)), ver_no(inn_dim2), ver_no_next(inn_dim2)) 

      cell_no = 0
      ver_no = 0
      d_val = 0


      count_T1 = 0
      do ic = 1, Lx*Ly
        vx = v(1, inn(1:num(ic), ic))
        vy = v(2, inn(1:num(ic), ic))

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
          len_d = dsqrt(dx**2 + dy**2)

          if(len_d < min_d_T1)then
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

      !  cellNoT1 = cell_no(chosen_index)
      !  verNoT1 = ver_no(chosen_index)
      !  verNoNextT1 = ver_no_next(chosen_index)
        
        cellNoT1 = cell_no(1)
        verNoT1 = ver_no(1)
        verNoNextT1 = ver_no_next(1)

      else
        cellNoT1 = 0
        verNoT1 = 0
        verNoNextT1 = 0
      end if

!      write(*,*)cellNoT1, verNoT1
      
   end subroutine find_T1_2

   subroutine find_T1_Affected_2
     implicit none

     call find_T1_2
    ! write(*,*)inn(verNoT1,cellNoT1), inn(verNoNextT1, cellNoT1)

     do ic = 1, Lx*Ly
       do jc = 1, num(ic)
         if((inn(jc,ic)==inn(verNoT1,cellNoT1)).or. & 
           inn(jc,ic)==inn(verNoNextT1,cellNoT1))then

           write(*,*)ic, jc

         end if

       end do

     end do



   end subroutine find_T1_Affected_2



end module T1_swap_new
