module allocation

   implicit none
!   save

      integer :: Nc
      integer :: nrun, nrun2_initialTime
      real*8 :: nrun2_initialTime_r
      integer*4 :: i,j,k,it,jp,jm,n,nn,totT
      integer*4 :: ii,jj, iki, ik, jk, kke
      integer :: ic, jc
      real*8 :: totTr
      real*8 :: min_d_T1, min_area_T2
      integer :: count_T1, cellNoT1, verNoT1, verNoNextT1
      integer :: count_T2, cellNoT2, verNoT2
      real*8, allocatable, dimension(:) :: Total_T1_count, Total_T2_count

      integer :: next_idx, prev_idx, chosen_index
      real*8 :: len_d,len_d_sq, dx, dy, rr
      integer*4, dimension(:), allocatable :: cell_no, ver_no, ver_no_next
      real*8, dimension(:), allocatable :: d_val
      real*8, dimension(:), allocatable :: area_val


      logical :: if_Do_T1, if_Do_T2

      real*8 :: len_d_T1
      integer*4 :: Lx,Ly

      integer*4, dimension(:), allocatable :: num
      real*8, dimension(:,:), allocatable:: v
      integer*4, dimension(:,:), allocatable :: inn
      integer*4, dimension(:), allocatable :: borderver
      integer*4 :: borderdim
      integer*4, allocatable, dimension(:) :: coordNum, bound
      integer :: boundary_count
   
   
      real*8 :: eta, lambda, beta, gamm, Ao, Co,  dt

      real*8, allocatable, dimension(:) :: coefficients
      real*8, dimension(:),  allocatable ::  fxx_temp, fyy_temp
      real*8, dimension(:), allocatable:: fxx, fyy
      real*8, dimension(:), allocatable:: fxx_ran, fyy_ran
      integer*4, dimension(:), allocatable :: workingzone
      real*8, dimension(:), allocatable :: vx,vy

      real*8, dimension(:), allocatable :: edgelength, edgelengthIn
      real*8 :: area, perimeter

      real*8 :: totalarea
      real*8 ::  pi, rann

      integer*4 :: num_dim, inn_dim1,inn_dim2, v_dim1, v_dim2

      integer*4, dimension(:), allocatable :: mainarea,leftP,rightP,topP,bottomP,corners,boundary, boundary_occ
      integer, dimension(:), allocatable :: inside1, inside2
      integer*4, dimension(:), allocatable :: GleftP,GrightP,GtopP,GbottomP,Gcorners,Gboundary

      real*8, allocatable, dimension(:) :: Energy, ShearStress
      character(100) :: fname_inn, fname_num, fname_v, fname_force
      character(100) :: fname_inn2, fname_num2, fname_v2
      integer*4 :: iunit_inn, iunit_num, iunit_v

      integer*4 :: it_dump, T1_time_interval, T2_time_interval

      real*8, allocatable, dimension(:) ::  mot, mot0
      real*8 :: etas_max, etas_min, mot_Lc
      real*8, allocatable, dimension(:) :: msdt, cellCentInit
      logical :: if_Shear_tissue
      real*8 :: sudden_shearStrength
      integer :: sudden_shearWhen, Oscl_shearWhen
      real*8 :: strainRate, Oscl_freq_wo, Oscl_shearStrength

      integer :: sin_perturb_when
      real*8 :: sin_perturb_strength, sin_perturb_waveNumber
      logical :: if_Perturb_tissue, if_sin_perturb, if_dirac_comb
      integer :: comb_onPeriod, comb_offPeriod


      real*8 :: sigma(2,2)
      real*8, allocatable, dimension(:,:) :: TotalSigma


      logical :: if_motility_gradient, if_Fixed_boundary, if_Sudden_Shearing, if_Oscillatory_Shearing
      logical :: if_motility_decay
      real*8 :: motility_decay_timeScale

      logical :: if_bottom_borders_fixed

      real*8, allocatable, dimension(:,:) :: cellcen
      real*8 :: global_cellCenX, global_cellCenY

      integer, allocatable, dimension(:)  :: bottom_border, top_border, left_border, right_border, all_borders
      integer :: bottom_border_count, top_border_count, left_border_count, right_border_count
      integer :: all_border_count
      
      integer, allocatable, dimension(:) :: vertex_occurance_count

      logical :: if_limb_force
      real*8 :: limb_force_strength

      logical :: if_motility_hotspot
      integer :: number_of_hotspot
      real*8 :: sigma_hotspot

      logical :: if_cell_division
      real*8 :: area_0

      integer, allocatable, dimension(:) :: chosen_cell
      integer :: chosen_cell_count


   contains 

   subroutine read_input

     implicit none
  
!     open(unit=121, file='para.in', status='old'); 
     open(112, file='para1_in.dat', status='old') 
     open(unit=121, file='para2_in.dat', status='old'); 
  
     read(121,*) Lx
     read(121,*) Ly
     read(121,*) num_dim
     read(121,*) v_dim1
     read(121,*) v_dim2
     read(121,*) inn_dim1
     read(121,*) inn_dim2

     read(112,*) nrun
     read(112,*) nrun2_initialTime_r
     read(112,*) Ao
     read(112,*) Co
     read(112,*) lambda
     read(112,*) beta
     read(112,*) gamm
     read(112,*) eta
     read(112,*) totTr
     read(112,*) dt
     read(112,*) if_Do_T1
     read(112,*) min_d_T1
     read(112,*) if_Do_T2
     read(112,*) min_area_T2
     read(112,*) if_Fixed_boundary
     read(112,*) if_bottom_borders_fixed
     read(112,*) it_dump
     read(112,*) T1_time_interval
     read(112,*) T2_time_interval
     read(112,*) if_motility_gradient
     read(112,*) etas_max
     read(112,*) etas_min
     read(112,*) mot_Lc
     read(112,*) if_motility_decay
     read(112,*) motility_decay_timeScale
     read(112,*) if_motility_hotspot
     read(112,*) number_of_hotspot
     read(112,*) sigma_hotspot
     read(112,*) if_Shear_tissue
     read(112,*) if_Sudden_Shearing
     read(112,*) sudden_shearStrength
     read(112,*) sudden_shearWhen
     read(112,*) if_Oscillatory_Shearing
     read(112,*) Oscl_shearStrength
     read(112,*) Oscl_shearWhen
     read(112,*) Oscl_freq_wo
     read(112,*) if_Perturb_tissue
     read(112,*) if_sin_perturb
     read(112,*) sin_perturb_when
     read(112,*) sin_perturb_strength
     read(112,*) sin_perturb_waveNumber
     read(112,*) if_dirac_comb
     read(112,*) comb_onPeriod
     read(112,*) comb_offPeriod
     read(112,*) if_limb_force
     read(112,*) limb_force_strength
     read(112,*) if_cell_division
     read(112,*) area_0

     




     totT = int(totTr)
     nrun2_initialTime = int(nrun2_initialTime_r)

     pi = acos(-1.0d0)

     close(121)
     close(112)
  
   end subroutine read_input

   subroutine allocate_arrays
     implicit none
     call read_input

     allocate(num(num_dim))
     allocate(v(v_dim1,v_dim2))
     allocate(inn(inn_dim1,inn_dim2))
     allocate(fxx(v_dim2), fyy(v_dim2))
     allocate(fxx_ran(v_dim2), fyy_ran(v_dim2))
     allocate(fxx_temp(v_dim2), fyy_temp(v_dim2))
     allocate(edgelength(Lx*Ly*inn_dim1))
     allocate(edgelengthIn(Lx*Ly*inn_dim1))
     allocate(workingzone(Lx*Ly+2*Lx+2*Ly+4))
     allocate(mainarea(Lx*Ly),leftP(Ly),rightP(Ly),topP(Lx),bottomP(Lx),corners(4))
     allocate(boundary(2*Lx + 2*Ly))
     allocate(GleftP(Ly),GrightP(Ly),GtopP(Lx),GbottomP(Lx),Gcorners(4))
     allocate(Gboundary(2*Lx + 2*Ly + 4))
     allocate(inside1((Lx - 2) * (Ly - 2)))
     allocate(inside2((Lx - 4) * (Ly - 4)))
     allocate(coefficients(6))
     allocate(Energy(totT))
     allocate(msdt(totT))
     allocate(Total_T1_count(totT), Total_T2_count(totT))
     allocate(cellCentInit(Lx*Ly-4*Lx - 4*Ly +16))
     allocate(borderver(v_dim2)) 
     allocate(mot(v_dim2),mot0(v_dim2))
     allocate(coordNum(Lx*Ly))
     allocate(bound(2*Lx + 2*Ly))
     allocate(cellcen(Lx*Ly, 2))

     allocate(d_val(Lx*Ly*6))
     allocate(area_val(Lx*Ly))
     allocate(ver_no(inn_dim2), ver_no_next(inn_dim2))
     allocate(bottom_border(Lx*6), top_border(Lx*6), left_border(Lx*6), right_border(Lx*6))
     allocate(all_borders(Lx*Ly*6))
     allocate(TotalSigma(2,2))
     allocate(ShearStress(totT))

     allocate(chosen_cell(v_dim2))

    end subroutine allocate_arrays

    subroutine read_data

      if(nrun.eq.1)then
       open(1, file='num_in.dat',status='old')
       open(2, file='v_in.dat', status='old')
       open(3, file='inn_in.dat', status='old')

       
       
    
       read(1,*)num
       read(2,*)v
       read(3,*)inn

    
    
       close(1)
       close(2)
       close(3)

     elseif(nrun.eq.2)then


       
        iunit_inn = 533
        iunit_num = 623
        iunit_v = 959

        write(fname_inn2, '("data/inn_", I8.8,".dat")')(nrun2_initialTime)
        write(fname_v2, '("data/v_", I8.8,".dat")')(nrun2_initialTime)
        write(fname_num2, '("data/num_", I8.8,".dat")')(nrun2_initialTime)
        write(*,*)fname_inn2

        open(unit = iunit_inn, file=fname_inn2, form = 'unformatted', status='old')
        open(unit = iunit_v, file=fname_v2, form = 'unformatted', status='old')
        open(unit = iunit_num, file=fname_num2, form = 'unformatted', status='old')


        read(iunit_inn)inn
        read(iunit_v)v
        read(iunit_num)num
       
        close(iunit_inn)
        close(iunit_num)
        close(iunit_v)





     end if

     coefficients = [lambda, beta, gamm, Ao, Co, eta]

     allocate(cell_no(maxval(inn)))   ! allocating here to increase speed
     allocate(vertex_occurance_count(maxval(inn)))


    end subroutine read_data

    subroutine write_output

      integer :: iunit_inn, iunit_num, iunit_v, iunit_force


       iunit_inn = 532
       iunit_num = 621
       iunit_v = 958
       iunit_force = 961
 
       if(nrun.eq.1)then
         write(fname_inn, '("data/inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/num_", I8.8,".dat")')(it)
         write(fname_v, '("data/v_", I8.8,".dat")')(it)
         write(fname_force, '("data/force_", I8.8,".dat")')(it)
       elseif(nrun.eq.2)then
         write(fname_inn, '("data/nrun2_inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/nrun2_num_", I8.8,".dat")')(it)
         write(fname_v, '("data/nrun2_v_", I8.8,".dat")')(it)
         write(fname_force, '("data/nrun2_force_", I8.8,".dat")')(it)
       end if



       open(unit = iunit_inn, file=fname_inn, form = 'unformatted', status='unknown')             
       open(unit = iunit_num,file=fname_num, form ='unformatted', status='unknown')
       open(unit = iunit_v, file=fname_v, form='unformatted', status='unknown')
       open(unit = iunit_force, file=fname_force, form='unformatted', status='unknown')
       ! open(unit = iunit_force, file=fname_force, status='unknown')
 
 
       write(iunit_inn)((inn(i,j),i=1,inn_dim1),j=1,inn_dim2)
       write(iunit_num)(num(i), i=1,num_dim)
       write(iunit_v)((v(i,j), i=1,v_dim1),j=1,v_dim2)
       write(iunit_force)(fxx(i), fyy(i), i = 1, v_dim2)
      
      
 
       if(it.eq.totT)then
         open(unit = 711, file='data/Energy.dat', form='unformatted',  status='unknown')
         write(711)Energy
         close(711)
       end if
       if(it.eq.totT)then
         open(unit = 715, file='data/ShearStress.dat', form='unformatted',  status='unknown')
         write(715)ShearStress
         close(715)
       end if
       if(it.eq.totT)then
         open(unit = 717, file='data/T1_count.dat', form='unformatted',  status='unknown')
         write(717)Total_T1_count
         close(717)
       end if

       if(it.eq.totT)then
         open(unit = 719, file='data/T2_count.dat', form='unformatted',  status='unknown')
         write(719)Total_T2_count
         close(719)
       end if


       if(it.eq.1)then
         open(unit=720, file='data/motility_store.dat', form='unformatted',status='unknown')
         write(720)(mot(i), i = 1, v_dim2)
         close(720)

       end if



 
 
       close(iunit_inn)
       close(iunit_num)
       close(iunit_v)
       close(iunit_force)

    end subroutine write_output





end module allocation
