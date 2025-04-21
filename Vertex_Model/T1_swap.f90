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
      do ic = 1, Lx*Ly

!        if(num(ic).le.3)cycle

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
     integer :: temp_array(Lx*Ly)
     

    ! write(*,*)inn(verNoT1,cellNoT1), inn(verNoNextT1, cellNoT1)
     

    ik = 0
     do ic = 1, Lx*Ly
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
     
    call PerpendicularBisector2(OldX,OldY,scaling_factor,NewX,NewY)

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
        
       
       if(verNo_in1_indx.ne.0.and.verNo_in2_indx.eq.0)then

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


       if(verNo_in2_indx.ne.0.and.verNo_in1_indx.eq.0)then
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
      
!       write(*,*)'inn_updated1',inn(:, cellNo)



     end if  ! Occurrence if
    
   end do   ! ik do
     



    end subroutine T1_core



    subroutine Do_T1


      implicit none
      integer :: il, im

      T1_pass = .true.
      
      call find_T1

      if(T1_pass)then
        if(if_Fixed_boundary)then
          call find_T1_Affected


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

      

      if(T1_pass.and.if_bottom_borders_fixed)then
      call Find_boundary_dynamic
        do im = 1, bottom_border_count
          if(bottom_border(im).eq.inn(verNoT1, cellNoT1) &
            .or.bottom_border(im).eq.inn(verNoT1+1, cellNoT1))then
            T1_pass = .false.
            write(*,*)'Bottom Border Ignored'
          end if
        end do
      end if


     
!      write(*,*)T1_pass,count_T1
      
      if(T1_pass)then
          call find_T1_Affected
          call T1_core    
          write(*,*)"T1 Happened at it = ", it
          Total_T1_count(it) = Total_T1_count(it) + 1
        end if



    end subroutine Do_T1






end module T1_swap
