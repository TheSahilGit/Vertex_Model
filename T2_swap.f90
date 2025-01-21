module T2_swap
  use array_info
  use Geometry
  use allocation



  logical :: T2_pass

  contains

    subroutine find_T2
      implicit none


      cell_no = 0
      ver_no = 0
      area_val = 0

      count_T2 = 0

      do ic = 1, Lx*Ly
        if(num(ic).le.3)then
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


!      if(count_T2>0)then
!        call random_number(rr)
!        chosen_index = int(rr*count_T2+1)
!        !cellNoT2 = cell_no(chosen_index)
!        cellNoT2 = 111  !cell_no(chosen_index)
!      else
!        T2_pass = .false.
!      end if


    end subroutine find_T2
    

    subroutine find_T2_Affected
      implicit none
      integer :: temp_array(Lx*Ly)

      integer :: ik

      !write(*,*)cellNoT2

!      ik = 0
!      do ic = 1,Lx*Ly
!        do jc = 1, num(ic)
!
!         ! write(*,*)inn(1,cellnoT2)
!          if(inn(jc,ic)==inn(1,cellnoT2).or. &
!            inn(jc,ic)==inn(2,cellnoT2) .or. &
!            inn(jc,ic)==inn(3,cellnoT2))then
!            ik = ik + 1
!            temp_array(ik) = ic
!          end if
!
!        end do
!      end do



    end subroutine find_T2_Affected
!
!    subroutine T2_core
!
!
!    end subroutine T2_core
!
!    subroutine do_T2
!
!    end subroutine do_T2




end module T2_swap
