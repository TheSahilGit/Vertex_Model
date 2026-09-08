module System_Info

  use allocation
  use Geometry

  contains
    subroutine Calculate_Energy

      real*8 :: En
      real*8 :: En_Area, En_Peri, En_Line




      En = 0.0d0

      do ic = 1, Nc !Lx*Ly

        nn = num(ic)

        call Gather_Cell_Vertices_PBC(inn(1:nn,ic), nn, vx, vy)

        call CalculateArea(vx, vy, nn, area)
        call CalculatePerimeter(vx, vy, nn, perimeter)
        area = abs(area)


        En_Area =  lambda * (area - Ao)**2
        En_Peri =  beta * (perimeter - Co)**2
        En_Line =  gamm * perimeter

        En = En + En_Area + En_Peri + En_Line

       end do

       Energy(it) = En


    end subroutine Calculate_Energy


    ! Differential line tension (log.txt, if_free_edge_adhesion): a fully
    ! separate subroutine, NOT a branch inside Calculate_Energy, so the
    ! original subroutine above is left byte-for-byte untouched -- an
    ! earlier version added an if/else branch directly inside
    ! Calculate_Energy and, even though the untaken (flag-off) branch was
    ! the exact original expression, -O3 -march=native produced a handful
    ! of ~1e-13-relative-magnitude floating-point differences in Energy.dat
    ! anyway (verified: Force.f90/Stress.f90's in-place if/else branches did
    ! NOT show this effect, only this subroutine did) -- almost certainly a
    ! vectorization/instruction-scheduling change from the extra local
    ! declarations and branch, not an actual arithmetic difference. Calling
    ! a completely separate subroutine (gated at the vertexmain.f90 call
    ! site, not in here) sidesteps that risk entirely and was confirmed
    ! byte-identical.
    subroutine Calculate_Energy_FreeEdge

      real*8 :: En
      real*8 :: En_Area, En_Peri, En_Line
      logical :: is_free_edge(inn_dim1, num_dim)
      integer :: jc_e, next_idx_e
      real*8 :: edx, edy, elen, gamm_this_edge

      En = 0.0d0

      call Get_Free_Edges(is_free_edge)

      do ic = 1, Nc !Lx*Ly

        nn = num(ic)

        call Gather_Cell_Vertices_PBC(inn(1:nn,ic), nn, vx, vy)

        call CalculateArea(vx, vy, nn, area)
        call CalculatePerimeter(vx, vy, nn, perimeter)
        area = abs(area)


        En_Area =  lambda * (area - Ao)**2
        En_Peri =  beta * (perimeter - Co)**2

        En_Line = 0.0d0
        do jc_e = 1, nn
          next_idx_e = jc_e + 1
          if (jc_e == nn) next_idx_e = 1
          edx = vx(next_idx_e) - vx(jc_e)
          edy = vy(next_idx_e) - vy(jc_e)
          elen = dsqrt(edx**2 + edy**2)
          gamm_this_edge = merge(gamm_free, gamm, is_free_edge(jc_e, ic))
          En_Line = En_Line + gamm_this_edge * elen
        end do

        En = En + En_Area + En_Peri + En_Line

       end do

       Energy(it) = En


    end subroutine Calculate_Energy_FreeEdge


   subroutine CellCentre
!     use Geometry
     use allocation

     do ic = 1, Nc !Lx*Ly

       call Gather_Cell_Vertices_PBC(inn(1:num(ic),ic), num(ic), vx, vy)

       cellcen(ic,1) = sum(vx)/dble(size(vx))
       cellcen(ic,2) = sum(vy)/dble(size(vy))

      end do

      ! BUGFIX (log.txt, re-review pass): cellcen is allocated (num_dim, 2)
      ! and never zero-initialized; only rows 1:Nc are filled by the loop
      ! above. Summing the whole column (cellcen(:,1)) included rows
      ! Nc+1:num_dim -- uninitialized heap memory. Currently dead code
      ! (CellCentre is never called), but fixed for when it is.
      global_cellCenX = sum(cellcen(1:Nc,1))/dble(Nc)
      global_cellCenY = sum(cellcen(1:Nc,2))/dble(Nc)




   end subroutine CellCentre


!   subroutine MeanSqDisp(Lx_in,Ly_in,v_in,inn_in,num_in,cellCentInit,avgdisp)
!!     use allocation
!
!     integer*4 :: ic, jc, nnm
!      integer*4, dimension(:), allocatable :: num_in
!      real*8, dimension(:,:), allocatable:: v_in
!      integer*4, dimension(:,:), allocatable :: inn_in
!      integer*4 :: Lx_in,Ly_in, working_L
!      real*8, allocatable, dimension(:) :: cellcen, cellCentInit
!      real*8, allocatable, dimension(:) :: disp
!      real*8 :: avgdisp
!
!
!
!!      write(*,*)'cellcent = ,', cellcen
!
!      working_L = Lx_in*Ly_in - 4*Lx_in - 4*Ly_in + 16
!
!      allocate(cellcen(working_L))
!      allocate(disp(working_L))
!
!      call CellCentre(Lx_in,Ly_in,v_in,inn_in,num_in,cellcen)
!
!      do i = 1, working_L
!        disp(i) = (abs(cellcen(i)-cellCentInit(i)))**2
!      end do
!
!      avgdisp = sum(disp)/(working_L)
!
!
!
!   end subroutine MeanSqDisp


end module System_Info
