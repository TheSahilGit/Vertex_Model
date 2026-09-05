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
      real*8 :: beta_0

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
      character(100) :: fname_Myosin, fname_cell_identity
      character(100) :: fname_inn2, fname_num2, fname_v2
      integer*4 :: iunit_inn, iunit_num, iunit_v
      ! T1/T2 event log (log.txt: spatial/cell-identity tracking for
      ! Movie_Code.m). Opened once (Open_Event_Logs) before the main loop,
      ! written to incrementally (one line per actual event, from Do_T1/
      ! Do_T2), closed once (Close_Event_Logs) after it.
      integer, parameter :: iunit_T1events = 941, iunit_T2events = 942

      integer*4 :: it_dump, T1_time_interval, T2_time_interval
      integer*4 :: summary_dump_interval

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

      logical :: if_motility

      logical :: if_motility_gradient, if_Fixed_boundary, if_Sudden_Shearing, if_Oscillatory_Shearing

      ! Eulerian (spatial-field) vs Lagrangian (per-vertex tag) motility
      ! gradient, per user request (log.txt): if_motility_gradient assigns
      ! mot once, at it=1, from each vertex's position at that instant
      ! (Give_Motility_Gradient, called only before the main loop) -- once
      ! assigned, a vertex carries that value wherever it physically moves
      ! for the rest of the run (Lagrangian). Under active motility forces,
      ! the highest-motility vertices (originally at the gradient's low end)
      ! also move the most, so the visible gradient erodes over the run
      ! (confirmed directly: corr(cell mean mot, cell's current y) dropped
      ! from -0.66 at it=20000 to -0.31 at it=200000 in a real run). Setting
      ! if_motility_Eulerian=.true. instead re-derives mot from each
      ! vertex's CURRENT position every step (see the Give_Motility_Gradient
      ! call in vertexmain.f90's main loop), so the gradient stays anchored
      ! to physical location -- a true spatial field -- instead of drifting
      ! with the tissue. Defaults .false. (Lagrangian) to keep every
      ! existing para_Simulation.dat's behavior unchanged; has no effect
      ! unless if_motility_gradient is also .true.
      logical :: if_motility_Eulerian

      logical :: if_motility_decay
      real*8 :: motility_decay_timeScale

      logical :: if_bottom_borders_fixed, if_top_borders_fixed

      ! True runtime periodic boundary conditions (log.txt, user request).
      ! Placed here, right after the other boundary-condition flags (its
      ! read(112,*) counterpart is inserted right after if_top_borders_
      ! fixed in read_input, mirrored at the same position in
      ! para_Simulation.dat) since if_PBC is fundamentally a fourth
      ! boundary-condition choice, mutually exclusive with the three
      ! above -- there is no free edge left for if_Fixed_boundary/
      ! if_bottom_borders_fixed/if_top_borders_fixed to pin once the
      ! tissue is periodic in both x and y (enforced in read_input).
      ! Requires a mesh built with if_periodic=.true. in
      ! para_MeshGen.dat/Generate_Initial_Mesh.f90 (a genuinely
      ! wrap-connected torus, not the default free-boundary mesh) --
      ! nothing here checks that at runtime; using if_PBC with an
      ! ordinary free-boundary mesh would silently treat its real free
      ! edge as if it wrapped, which is wrong.
      logical :: if_PBC
      ! Exact periodic box lengths, set once from Lx,Ly (cast to real*8)
      ! at the end of read_data, once Lx,Ly are known -- deliberately a
      ! separate name from the integer lattice-count Lx,Ly (used all over
      ! the rest of the code for index-range/array-sizing logic) to avoid
      ! any accidental mixing of "cell count" and "physical length"
      ! meanings. The mesh generator's ghost-copy tiling places periodic
      ! images at exactly +-Lx,+-Ly, so this is exact, not approximate.
      real*8 :: Lx_box, Ly_box

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

      logical :: if_squeeze_tissue
      integer :: squeeze_when
      real*8 :: percent_squeeze


      logical :: if_active_contractility
      real*8 :: tau_contr, active_contr_strength


      logical :: if_ABP
      real*8 :: vo, Dr
      real*8, dimension(:), allocatable :: fxx_ABP, fyy_ABP, theta_ABP, rot_noise


      logical :: if_RhoROCK, if_RK4, if_myosin_noise, if_gaussian_noise
      real*8, dimension(:), allocatable :: Rho, ROCK, Myosin
      real*8 :: A_Rho, nhill, K_hill, A_ROCK, A_Myosin, D_rho, D_ROCK, D_Myosin
      real*8 :: Myosin_Coupling_Strength, Myosin_noise_strength
      logical :: if_coupling_noise
      real*8 :: coupling_noise_strength

      character(500), dimension(:), allocatable :: cell_identity

      logical :: if_polar_motility
      real*8 :: polar_motility_strength
      real*8, dimension(:), allocatable :: fxx_Polar, fyy_Polar

      integer :: n_inside, n_outside
      real*8 :: radius_from_core
      real*8, allocatable, dimension(:,:) :: cell_centers
      integer, allocatable, dimension(:) :: inside_cells, outside_cells
      
        



   contains 

   subroutine read_input

     implicit none
  
!     open(unit=121, file='para.in', status='old'); 
     open(112, file='para_Simulation.dat', status='old')
     open(unit=121, file='para_MeshDims.dat', status='old');
  
     read(121,*) Lx
     read(121,*) Ly
     read(121,*) num_dim
     read(121,*) v_dim1
     read(121,*) v_dim2
     read(121,*) inn_dim1
     read(121,*) inn_dim2

     ! Exact periodic box lengths -- see the Lx_box/Ly_box declaration
     ! comment. Computed here (Lx,Ly are already known) rather than at the
     ! end of read_input, so it's available regardless of where if_PBC
     ! itself is later validated.
     Lx_box = dble(Lx)
     Ly_box = dble(Ly)

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
     read(112,*) if_top_borders_fixed
     read(112,*) if_PBC

     ! Validation (log.txt): if_PBC means the tissue is periodic in both x
     ! and y -- there is no free edge left, so a flag that means "pin/
     ! detect a free edge" is meaningless (and, worse, would silently
     ! mispin/misdetect on a periodic mesh's genuine interior vertices).
     ! Stop rather than run with an ambiguous/contradictory configuration.
     if (if_PBC) then
       if (if_Fixed_boundary .or. if_bottom_borders_fixed .or. if_top_borders_fixed) then
         write(*,*) 'read_input: if_PBC is incompatible with if_Fixed_boundary/', &
           'if_bottom_borders_fixed/if_top_borders_fixed (no free edge exists under PBC).'
         stop 1
       end if
     end if

     read(112,*) it_dump
     read(112,*) T1_time_interval
     read(112,*) T2_time_interval
     read(112,*) if_motility
     read(112,*) if_motility_gradient
     read(112,*) if_motility_Eulerian
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
     read(112,*) if_squeeze_tissue
     read(112,*) squeeze_when
     read(112,*) percent_squeeze
     read(112,*) if_active_contractility
     read(112,*) tau_contr
     read(112,*) active_contr_strength
     read(112,*) if_ABP
     read(112,*) vo
     read(112,*) Dr
     read(112,*) if_RhoROCK
     read(112,*) if_RK4
     read(112,*) nhill
     read(112,*) K_hill 
     read(112,*) A_Rho
     read(112,*) A_ROCK 
     read(112,*) A_Myosin
     read(112,*) D_rho
     read(112,*) D_ROCK
     read(112,*) D_Myosin
     read(112,*) Myosin_Coupling_Strength
     read(112,*) if_myosin_noise
     read(112,*) if_gaussian_noise
     read(112,*) Myosin_noise_strength
     read(112,*) if_coupling_noise
     read(112,*) coupling_noise_strength
     read(112,*) if_polar_motility
     read(112,*) polar_motility_strength

     



     beta_0 = beta

     ! BUGFIX (log.txt): beta/gamm must be nondimensionalized exactly ONCE here.
     ! Previously this rescale was done inside Force_Calculation (Force.f90) and
     ! again inside Calculate_StressTensor (Stress.f90), both of which run every
     ! timestep and mutate the global beta/gamm in place -- causing beta/gamm to
     ! be divided by (lambda*Ao) again on every call (exponential decay/blow-up
     ! over the run). Doing it once here, right after read, fixes that.
     beta = beta/(lambda*Ao)
     gamm = gamm/(lambda*(Ao)**1.5)

     totT = int(totTr)
     nrun2_initialTime = int(nrun2_initialTime_r)

     ! CHANGE (log.txt, MATLAB live-plotting request): Energy/ShearStress/
     ! T1_count/T2_count are now rewritten periodically during the run (not
     ! just once at the end) so they can be plotted while a simulation is
     ! still in progress -- see write_output. Naively rewriting the whole
     ! growing (1:it) slice every it_dump would cost O(totT^2/it_dump) total
     ! disk I/O over the run: for this session's actual production config
     ! (totT=1e7, it_dump=100) that works out to ~4TB written for Energy.dat
     ! alone (~16TB across all four files) -- large enough to dominate the
     ! run's wall-clock time. Instead, rewrite them only every
     ! summary_dump_interval, a multiple of it_dump chosen so there are
     ! about 200 rewrites total over the whole run regardless of totT/
     ! it_dump: this bounds the cumulative I/O to a small, fixed multiple of
     ! the final file size (~100x), while still refreshing often enough to
     ! watch a long run's progress. Never coarser than it_dump itself (a run
     ! with very few total dumps just gets a summary write every dump, same
     ! as before).
     summary_dump_interval = max(it_dump, it_dump * nint(dble(totT) / dble(it_dump * 200)))

     pi = acos(-1.0d0)

     ! Further if_PBC validation (log.txt) -- these flags are read later in
     ! this same subroutine than if_Fixed_boundary/if_bottom_borders_fixed/
     ! if_top_borders_fixed above, so they're checked here instead:
     !  - if_squeeze_tissue/if_limb_force both assume a bounded domain with
     !    a real bottom edge / real corners to act on -- meaningless under
     !    full periodicity (no edge, no corners).
     !  - if_motility_gradient's Eulerian mode (if_motility_Eulerian)
     !    recomputes mot = etas_max*exp(-y/(mot_Lc*Ly)) from each vertex's
     !    CURRENT y every step -- not periodic in y, so a vertex wrapping
     !    y=Ly->0 would see a discontinuous motility jump. Disallowed for
     !    now per explicit user decision; a periodic-safe gradient formula
     !    is a separate future feature, not attempted here.
     if (if_PBC) then
       if (if_squeeze_tissue) then
         write(*,*) 'read_input: if_PBC is incompatible with if_squeeze_tissue', &
           ' (assumes a real bottom edge to compress against).'
         stop 1
       end if
       if (if_limb_force) then
         write(*,*) 'read_input: if_PBC is incompatible with if_limb_force', &
           ' (assumes real tissue corners to push outward from).'
         stop 1
       end if
       if (if_motility_gradient .and. if_motility_Eulerian) then
         write(*,*) 'read_input: if_PBC with if_motility_gradient+if_motility_Eulerian', &
           ' is not supported yet (the gradient formula is not periodic in y).'
         stop 1
       end if
       if (if_Shear_tissue) then
         write(*,*) 'read_input: if_PBC with if_Shear_tissue is not supported yet -- ', &
           'ShearTissue directly displaces every vertex, which is not equivalent to ', &
           'genuine Lees-Edwards periodic shear (a separate, not-yet-implemented feature).'
         stop 1
       end if
     end if

     close(121)
     close(112)
  
   end subroutine read_input

   subroutine allocate_arrays
     implicit none
     integer :: iic
     call read_input

     allocate(num(num_dim))
     allocate(v(v_dim1,v_dim2))
     allocate(inn(inn_dim1,inn_dim2))
     allocate(fxx(v_dim2), fyy(v_dim2))
     allocate(fxx_ran(v_dim2), fyy_ran(v_dim2))
     allocate(fxx_Polar(v_dim2), fyy_Polar(v_dim2))
     allocate(fxx_temp(v_dim2), fyy_temp(v_dim2))
     allocate(fxx_ABP(v_dim2), fyy_ABP(v_dim2), theta_ABP(v_dim2), rot_noise(v_dim2))
     allocate(edgelength(Lx*Ly*inn_dim1))
     allocate(edgelengthIn(Lx*Ly*inn_dim1))
     allocate(workingzone(Lx*Ly+2*Lx+2*Ly+4))
     allocate(mainarea(Lx*Ly),leftP(Ly),rightP(Ly),topP(Lx),bottomP(Lx),corners(4))
     ! BUGFIX (log.txt, re-review pass): Get_Boundary_info actually builds
     ! boundary = [leftP, rightP, topP, bottomP, corners], i.e. 2*Lx+2*Ly+4
     ! elements (corners was missing from this initial sizing) -- harmless
     ! today only because "boundary" is allocatable and Fortran's automatic
     ! reallocation-on-assignment (gfortran default, no -fno-realloc-lhs in
     ! compile.sh) silently resizes it on that whole-array assignment.
     ! Sized correctly here so it isn't fragile if that assignment pattern
     ! ever changes to a partial/indexed one (as the border arrays already do).
     allocate(boundary(2*Lx + 2*Ly + 4))
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
     allocate(cellcen(num_dim, 2))

     ! BUGFIX (log.txt): d_val/ver_no/ver_no_next/area_val must all share one
     ! consistent worst-case capacity. Previously d_val was sized Lx*Ly*6 while
     ! ver_no/ver_no_next were sized only inn_dim2 (much smaller) -- if find_T1
     ! (T1_swap.f90) ever found more than inn_dim2 short edges in one pass,
     ! ver_no/ver_no_next overflowed. Also, sizing off the *initial* Lx*Ly (not
     ! num_dim, the reserved max-cell capacity) meant these silently became too
     ! small again once Nc grows via cell division. Now all sized from num_dim.
     ! BUGFIX (log.txt, re-review pass): the earlier num_dim*6 factor assumed
     ! an "average hexagonal cell" (~6 vertices/cell) but is not a proven
     ! bound against the actually configured per-cell vertex cap (inn_dim1);
     ! the true worst case (every live cell at inn_dim1 vertices, every edge
     ! short) is num_dim*inn_dim1. Sized off that instead.
     allocate(d_val(num_dim*inn_dim1))
     allocate(area_val(num_dim))
     allocate(ver_no(num_dim*inn_dim1), ver_no_next(num_dim*inn_dim1))
     ! BUGFIX (log.txt, re-review pass): left_border/right_border collect
     ! vertices along the VERTICAL (Ly) edges of the tissue, not the
     ! horizontal (Lx) ones -- they were sized Lx*6 like bottom/top_border,
     ! which only happened to be safe because the shipped mesh is square
     ! (Lx==Ly). Sized off Ly instead, matching the dimension they actually
     ! scale with (c.f. Get_Boundary_info's GleftP/GrightP, already Ly-sized).
     allocate(bottom_border(Lx*6), top_border(Lx*6), left_border(Ly*6), right_border(Ly*6))
     allocate(all_borders(Lx*Ly*6))
     allocate(TotalSigma(2,2))
     allocate(ShearStress(totT))
     ! BUGFIX (log.txt): ShearStress is unconditionally written to
     ! data/ShearStress.dat in write_output regardless of if_Shear_tissue.
     ! Without shearing it was only ever allocated, never assigned, so the
     ! output file contained uninitialized memory. Zero it here.
     ShearStress = 0.0d0

     allocate(chosen_cell(v_dim2))

     allocate(Rho(num_dim), ROCK(num_dim), Myosin(num_dim))

     allocate(cell_identity(num_dim))

     allocate(cell_centers(num_dim, 2), inside_cells(num_dim), outside_cells(num_dim))


     fxx = 0.0d0; fyy = 0.0d0
     fxx_ran = 0.0d0; fyy_ran = 0.0d0
     fxx_ABP = 0.0d0; fyy_ABP = 0.0d0; theta_ABP = 0.0d0
     fxx_Polar = 0.0d0; fyy_Polar = 0.0d0
     rot_noise = 0.0d0

     Rho = 0.0d0; ROCK = 0.0d0; Myosin = 0.0d0

      
     cell_identity = 'cell_0'
     do iic = 1, Lx*Ly
      write(cell_identity(iic), '(A,I0)') 'cell_', iic
     end do

     inside_cells = 0
     outside_cells = 0

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

     ! BUGFIX (log.txt, re-review pass): cell_no and vertex_occurance_count
     ! were sized off maxval(inn) -- the largest vertex index present in the
     ! *initial* mesh only, since read_data runs once before the main loop.
     ! vertex_occurance_count is indexed DIRECTLY by vertex id
     ! (Find_boundary_dynamic, Geometry.f90), and cell_division legitimately
     ! creates vertices with indices beyond that initial snapshot (up to
     ! v_dim2, the real reserved capacity) -- confirmed live: this produced
     ! "Index '4568' ... above upper bound of 4564" once cell division plus
     ! a T2 event grew the mesh past the initial maxval(inn). cell_no shares
     ! the same fragile sizing convention (it's filled in lockstep with
     ! ver_no in find_T1/find_T2's candidate scan) even though its true
     ! required capacity is candidate-count-based, not vertex-id-based.
     ! Sized off v_dim2 (the true reserved vertex capacity) for both.
     allocate(cell_no(v_dim2))
     allocate(vertex_occurance_count(v_dim2))

     ! BUGFIX (log.txt): Nc (live cell count) used to be hardcoded to Lx*Ly
     ! in vertexmain.f90 regardless of nrun. That's correct for nrun=1 (a
     ! freshly generated mesh always has exactly Lx*Ly live cells), but
     ! wrong for nrun=2: cell division (Proliferation.f90: Nc=Nc+1) and T2
     ! extrusion (T2_swap.f90: Nc=Nc-1) both change the true live-cell
     ! count over the course of a run, and every physics loop (Force.f90,
     ! T1_swap.f90, T2_swap.f90) is bounded by `do ic=1,Nc`. Restarting
     ! (nrun=2) from a snapshot taken *after* any division/T2 event had a
     ! true Nc different from Lx*Ly, yet Nc was silently reset to Lx*Ly on
     ! restart -- confirmed live: restarting from a snapshot with 608 live
     ! cells (576 initial + 32 divisions) came back reporting Nc=576, so
     ! every already-divided daughter cell would have been silently
     ! ignored by physics from that point on despite loading correctly
     ! into v/inn/num. Derived here instead, from the invariant (already
     ! relied on by T2_swap.f90 and every MATLAB analysis script) that
     ! cells 1:Nc are live and contiguous and num(k)=0 for every k>Nc --
     ! this gives the same Lx*Ly for nrun=1 (all of num(1:Lx*Ly) are
     ! nonzero, the rest zero, by construction) and the correct restored
     ! value for nrun=2.
     Nc = count(num /= 0)

     ! PBC sanity check (log.txt -- confirmed live by the user hitting
     ! exactly this): a mesh built with if_periodic=.true.
     ! (Generate_Initial_Mesh.f90) has every vertex touched by exactly 3
     ! cells -- zero free-boundary (degree <3) vertices, by construction
     ! (the same invariant that generator itself verifies before ever
     ! writing the mesh out). Running such a mesh with if_PBC=.false. feeds
     ! cells that genuinely straddle the periodic wrap into non-wrap-aware
     ! geometry (CalculateArea/Force_Calculation/etc treat their raw,
     ! far-apart stored coordinates as one ordinary polygon) -- this
     ! reliably corrupts area/perimeter for those cells into huge, wrong
     ! values, blowing up forces into NaN within a handful of steps.
     ! Detect the mismatch here, before the run starts, instead of letting
     ! it silently NaN partway through.
     block
       integer, allocatable :: deg_check(:)
       integer :: kk, ll
       logical :: mesh_looks_periodic

       allocate(deg_check(v_dim2))
       deg_check = 0
       do kk = 1, Nc
         do ll = 1, num(kk)
           deg_check(inn(ll,kk)) = deg_check(inn(ll,kk)) + 1
         end do
       end do

       mesh_looks_periodic = .true.
       do kk = 1, v_dim2
         if (deg_check(kk) > 0 .and. deg_check(kk) < 3) then
           mesh_looks_periodic = .false.
           exit
         end if
       end do
       deallocate(deg_check)

       if ((.not. if_PBC) .and. mesh_looks_periodic) then
         write(*,*) 'read_data: STOPPING -- this mesh has ZERO free-boundary vertices', &
           ' (every vertex has degree exactly 3) -- it looks like a periodic mesh', &
           ' (if_periodic=.true. in Generate_Initial_Mesh.f90) but if_PBC is .false. in', &
           ' para_Simulation.dat. Running non-PBC-aware physics on a wrap-connected mesh', &
           ' corrupts area/perimeter for any cell straddling the wrap and reliably blows', &
           ' up into NaN within a few steps. Set if_PBC=.true. to match this mesh, or', &
           ' generate a non-periodic mesh (if_periodic=.false.) to match if_PBC=.false.'
         stop 1
       end if

       if (if_PBC .and. (.not. mesh_looks_periodic)) then
         write(*,*) 'read_data: WARNING -- if_PBC is .true. but this mesh has free-boundary', &
           ' (degree < 3) vertices -- it does not look like a mesh built with', &
           ' if_periodic=.true. in Generate_Initial_Mesh.f90. This is not immediately', &
           ' dangerous (if_Fixed_boundary is already forced off under if_PBC, and no real', &
           ' vertex pair is far enough apart to trigger a spurious wrap), but the tissue''s', &
           ' true free edge is now just an ordinary unpinned edge, not a periodic one --', &
           ' if_PBC is not doing what you likely intend with this mesh.'
       end if
     end block

    end subroutine read_data

    subroutine Open_Event_Logs
      ! T1/T2 event log (log.txt): opened ONCE before the main it-loop
      ! (unlike write_output's snapshot files, which are reopened every
      ! dump) -- events are sparse and irregular, so they're appended to
      ! these two long-lived file handles as they happen, closed once at
      ! the very end (Close_Event_Logs). Same nrun2_ prefix convention as
      ! every other output file: a restart run's log never overwrites the
      ! original run's.
      implicit none
      character(100) :: fname_T1events, fname_T2events

      if (nrun.eq.1) then
        fname_T1events = 'data/T1_events.dat'
        fname_T2events = 'data/T2_events.dat'
      else
        fname_T1events = 'data/nrun2_T1_events.dat'
        fname_T2events = 'data/nrun2_T2_events.dat'
      end if

      ! BINARY (log.txt): one Fortran unformatted record per event (each
      ! WRITE statement on Do_T1/Do_T2's already-open unit produces its own
      ! fixed-size record, leading/trailing 4-byte length markers included
      ! -- the same convention every other output file in this codebase
      ! uses, just one small record per event instead of one big record for
      ! the whole run). All-real*8, matching this codebase's own habit of
      ! using real*8 even for integer-valued quantities (e.g.
      ! Total_T1_count) -- avoids ever mixing types within one unformatted
      ! record. Cell identities are stored as the numeric suffix of their
      ! 'cell_<N>' string (see CellIdNum below), 0.0d0 for 'cell_0'/no-cell
      ! padding -- matching this codebase's existing "0 means empty" idiom
      ! (e.g. num(i)==0). Roughly half the size of an equivalent formatted
      ! (ASCII) log, and keeps the same "readable while still running"
      ! property every other incrementally-written file here has.
      open(unit=iunit_T1events, file=fname_T1events, form='unformatted', status='replace')
      open(unit=iunit_T2events, file=fname_T2events, form='unformatted', status='replace')

    end subroutine Open_Event_Logs

    function CellIdNum(id_str) result(num_out)
      ! T1/T2 event log (log.txt): cell_identity is always exactly
      ! 'cell_'//<integer> (allocation.f90's read_data, T2_core's 'cell_0'
      ! placeholder) -- extract that integer directly rather than storing
      ! the whole padded string, for the compact binary event log above.
      implicit none
      character(len=*), intent(in) :: id_str
      real*8 :: num_out
      integer :: n, ios
      read(id_str(6:), *, iostat=ios) n
      if (ios /= 0) then
        num_out = 0.0d0
      else
        num_out = dble(n)
      end if
    end function CellIdNum

    subroutine Close_Event_Logs
      implicit none
      close(iunit_T1events)
      close(iunit_T2events)
    end subroutine Close_Event_Logs

    subroutine write_output

      integer :: iunit_inn, iunit_num, iunit_v, iunit_force
      integer :: iunit_Myosin
      integer :: iunit_cell_identity
      character(100) :: fname_Energy, fname_ShearStress, fname_T1count, &
        fname_T2count, fname_motility


       iunit_inn = 532
       iunit_num = 621
       iunit_v = 958
       iunit_force = 961
       iunit_Myosin = 966
       iunit_cell_identity = 967

       if(nrun.eq.1)then
         write(fname_inn, '("data/inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/num_", I8.8,".dat")')(it)
         write(fname_v, '("data/v_", I8.8,".dat")')(it)
         write(fname_force, '("data/force_", I8.8,".dat")')(it)
         write(fname_Myosin, '("data/Myosin_", I8.8,".dat")')(it)
         write(fname_cell_identity, '("data/cell_identity_", I8.8,".dat")')(it)
         fname_Energy = 'data/Energy.dat'
         fname_ShearStress = 'data/ShearStress.dat'
         fname_T1count = 'data/T1_count.dat'
         fname_T2count = 'data/T2_count.dat'
         fname_motility = 'data/motility_store.dat'
       elseif(nrun.eq.2)then
         write(fname_inn, '("data/nrun2_inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/nrun2_num_", I8.8,".dat")')(it)
         write(fname_v, '("data/nrun2_v_", I8.8,".dat")')(it)
         write(fname_force, '("data/nrun2_force_", I8.8,".dat")')(it)
         write(fname_Myosin, '("data/nrun2_Myosin_", I8.8,".dat")')(it)
         write(fname_cell_identity, '("data/nrun2_cell_identity_", I8.8,".dat")')(it)
         ! BUGFIX (log.txt): these 5 were still hardcoded to the same
         ! filenames as nrun=1 regardless of nrun -- unlike every other
         ! output above, an nrun=2 restart run silently overwrote the
         ! original nrun=1 run's Energy/ShearStress/T1_count/T2_count/
         ! motility_store files with its own (shorter, restarted) data,
         ! destroying that part of the original run's record. Matched to
         ! the same nrun2_ prefix convention as everything else here.
         fname_Energy = 'data/nrun2_Energy.dat'
         fname_ShearStress = 'data/nrun2_ShearStress.dat'
         fname_T1count = 'data/nrun2_T1_count.dat'
         fname_T2count = 'data/nrun2_T2_count.dat'
         fname_motility = 'data/nrun2_motility_store.dat'
       end if



       open(unit = iunit_inn, file=fname_inn, form = 'unformatted', status='unknown')             
       open(unit = iunit_num,file=fname_num, form ='unformatted', status='unknown')
       open(unit = iunit_v, file=fname_v, form='unformatted', status='unknown')
       open(unit = iunit_force, file=fname_force, form='unformatted', status='unknown')
       open(unit = iunit_Myosin,file=fname_Myosin, form ='unformatted', status='unknown')
       ! open(unit = iunit_force, file=fname_force, status='unknown')
       open(unit = iunit_cell_identity,file=fname_cell_identity, form ='unformatted', status='unknown')
 
 
       write(iunit_inn)((inn(i,j),i=1,inn_dim1),j=1,inn_dim2)
       write(iunit_num)(num(i), i=1,num_dim)
       write(iunit_v)((v(i,j), i=1,v_dim1),j=1,v_dim2)
       write(iunit_force)(fxx(i), fyy(i), fxx_ran(i), fyy_ran(i), & 
        fxx_ABP(i), fyy_ABP(i), & 
        fxx_Polar(i), fyy_Polar(i), i = 1, v_dim2)
       write(iunit_Myosin)(Rho(i), ROCK(i), Myosin(i), i = 1, num_dim)
       write(iunit_cell_identity)(cell_identity(i),  i=1,num_dim)
      
      
 
       ! CHANGE (log.txt, MATLAB live-plotting request): these four used to
       ! only be written once, at it==totT -- meaning Energy/ShearStress/
       ! T1_count/T2_count did not exist at all until a run fully finished,
       ! which blocked plotting them while a run is still in progress. Now
       ! rewritten periodically (every summary_dump_interval, computed in
       ! read_input -- capped to ~200 rewrites over the whole run regardless
       ! of totT/it_dump, to avoid an O(totT^2/it_dump) I/O blowup; see the
       ! comment at summary_dump_interval's computation) and always on the
       ! final timestep, each time writing only the (1:it) slice actually
       ! computed so far (entries it+1:totT are not yet valid). The file
       ! simply grows longer at each such write; MATLAB's
       ! fread(fid,Inf,'float64') pattern (LoadData.m and friends) already
       ! reads "however much is there" with no hardcoded length assumption,
       ! so this needs no MATLAB-side change.
       if (modulo(it, summary_dump_interval).eq.0 .or. it.eq.totT) then

         open(unit = 711, file=fname_Energy, form='unformatted',  status='unknown')
         write(711)Energy(1:it)
         close(711)

         open(unit = 715, file=fname_ShearStress, form='unformatted',  status='unknown')
         write(715)ShearStress(1:it)
         close(715)

         open(unit = 717, file=fname_T1count, form='unformatted',  status='unknown')
         write(717)Total_T1_count(1:it)
         close(717)

         open(unit = 719, file=fname_T2count, form='unformatted',  status='unknown')
         write(719)Total_T2_count(1:it)
         close(719)

       end if


       if(it.eq.1)then
         open(unit=720, file=fname_motility, form='unformatted',status='unknown')
         write(720)(mot(i), i = 1, v_dim2)
         close(720)

       end if



 
 
       close(iunit_inn)
       close(iunit_num)
       close(iunit_v)
       close(iunit_force)
       ! BUGFIX (log.txt): these two were never closed at all -- unlike
       ! inn/num/v/force just above, iunit_Myosin/iunit_cell_identity were
       ! left open after every single write_output call, relying on the
       ! NEXT call's open() (same fixed unit number, new filename) to
       ! implicitly close them. That leaves each snapshot's Myosin/
       ! cell_identity file in a not-yet-finalized state for up to a full
       ! it_dump interval (until the next dump event reopens the unit) --
       ! a real, independent bug from the read-side race this was found
       ! alongside (see that entry for the MATLAB-side half of the fix).
       close(iunit_Myosin)
       close(iunit_cell_identity)

    end subroutine write_output





end module allocation
