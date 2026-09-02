! Generate_Initial_Mesh.f90
!
! Pure-Fortran replacement for Main.m's MATLAB-based initial-mesh generation
! (see log.txt for the full writeup). Produces v_in.dat / inn_in.dat /
! num_in.dat / para_MeshDims.dat in exactly the format allocation.f90's
! read_input/read_data already expect -- no other Fortran file changes.
!
! Algorithm, mirroring Main.m's MeshGenerator/StoreData:
!   1. Hexagonal-offset (brick) lattice of Lx*Ly seed points.
!   2. Anisotropic random jitter: a single random draw per point, reused
!      (at different scales) for both its x and y offset -- this reproduces
!      Main.m's `ran` reuse exactly, including the resulting x/y-correlated
!      jitter direction, rather than "fixing" it into independent draws.
!   3. 8 periodic ghost copies (+-Lx, +-Ly and combinations) added purely so
!      edge cells clip against real neighbors instead of being unbounded --
!      this is a boundary-conditioning trick, not true periodic BCs.
!   4. Each real cell's Voronoi polygon is built directly via incremental
!      half-plane clipping against nearby candidate points (real + ghost),
!      processed in increasing distance order with an exact early-
!      termination test -- no general Delaunay/Voronoi algorithm needed,
!      since the point set is a mildly-perturbed lattice.
!   5. Vertices shared between cells are merged into one global list via a
!      uniform spatial grid (bucket-and-match within a small epsilon).
!   6. Each cell's vertex list is rotated to start at its min-x vertex, then
!      forced clockwise via a shoelace-sign check -- identical convention to
!      StoreData.m (Force.f90's area-gradient formula assumes clockwise
!      winding).
!
! Compiled standalone (see compile.sh) -- does not touch vertexmain.exe.

module spatial_grid
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   type grid_t
      integer :: nx, ny
      real(dp) :: xmin, ymin, cell
      integer, allocatable :: head(:,:)
      integer, allocatable :: nxt(:)
   end type grid_t

contains

   subroutine grid_init(g, xmin, ymin, xmax, ymax, cellsize, maxpts)
      type(grid_t), intent(out) :: g
      real(dp), intent(in) :: xmin, ymin, xmax, ymax, cellsize
      integer, intent(in) :: maxpts
      g%cell = cellsize
      g%xmin = xmin
      g%ymin = ymin
      g%nx = max(1, int((xmax - xmin) / cellsize) + 2)
      g%ny = max(1, int((ymax - ymin) / cellsize) + 2)
      allocate(g%head(g%nx, g%ny))
      g%head = 0
      allocate(g%nxt(maxpts))
      g%nxt = 0
   end subroutine grid_init

   subroutine grid_bucket(g, x, y, ix, iy)
      type(grid_t), intent(in) :: g
      real(dp), intent(in) :: x, y
      integer, intent(out) :: ix, iy
      ix = min(max(int((x - g%xmin) / g%cell) + 1, 1), g%nx)
      iy = min(max(int((y - g%ymin) / g%cell) + 1, 1), g%ny)
   end subroutine grid_bucket

   subroutine grid_insert(g, x, y, idx)
      type(grid_t), intent(inout) :: g
      real(dp), intent(in) :: x, y
      integer, intent(in) :: idx
      integer :: ix, iy
      call grid_bucket(g, x, y, ix, iy)
      g%nxt(idx) = g%head(ix, iy)
      g%head(ix, iy) = idx
   end subroutine grid_insert

   ! Returns candidate point indices whose bucket lies within `rad` of
   ! (x,y) -- caller must still check true distance <= rad; this only
   ! narrows the search.
   subroutine grid_query(g, x, y, rad, out_idx, out_n, maxout)
      type(grid_t), intent(in) :: g
      real(dp), intent(in) :: x, y, rad
      integer, intent(out) :: out_idx(:)
      integer, intent(out) :: out_n
      integer, intent(in) :: maxout
      integer :: ix0, ix1, iy0, iy1, ix, iy, p
      ix0 = min(max(int((x - rad - g%xmin) / g%cell) + 1, 1), g%nx)
      ix1 = min(max(int((x + rad - g%xmin) / g%cell) + 1, 1), g%nx)
      iy0 = min(max(int((y - rad - g%ymin) / g%cell) + 1, 1), g%ny)
      iy1 = min(max(int((y + rad - g%ymin) / g%cell) + 1, 1), g%ny)
      out_n = 0
      do iy = iy0, iy1
         do ix = ix0, ix1
            p = g%head(ix, iy)
            do while (p /= 0)
               out_n = out_n + 1
               if (out_n > maxout) stop 'grid_query: increase maxout'
               out_idx(out_n) = p
               p = g%nxt(p)
            end do
         end do
      end do
   end subroutine grid_query

end module spatial_grid


module mesh_gen_utils
   use spatial_grid
   implicit none
   integer, parameter :: MAXV = 64   ! generous per-cell polygon vertex cap

contains

   ! Clip convex polygon (px,py)(1:n), given in P-CENTERED coordinates (i.e.
   ! the cell's seed point P is the local origin), by the half-plane of
   ! candidate point q (also P-centered): keep points closer to the origin
   ! than to q. Standard Sutherland-Hodgman single-half-plane clip.
   subroutine clip_by_point(px, py, n, qx, qy)
      real(dp), intent(inout) :: px(MAXV), py(MAXV)
      integer, intent(inout) :: n
      real(dp), intent(in) :: qx, qy
      real(dp) :: qq, fi, fip1, t, newx(MAXV), newy(MAXV)
      integer :: i, ip1, m
      logical :: insidei, insideip1

      qq = 0.5d0 * (qx * qx + qy * qy)
      m = 0
      do i = 1, n
         ip1 = mod(i, n) + 1
         fi = px(i) * qx + py(i) * qy - qq
         fip1 = px(ip1) * qx + py(ip1) * qy - qq
         insidei = (fi <= 0.0d0)
         insideip1 = (fip1 <= 0.0d0)
         if (insidei) then
            m = m + 1
            if (m > MAXV) stop 'clip_by_point: MAXV exceeded'
            newx(m) = px(i)
            newy(m) = py(i)
         end if
         if (insidei .neqv. insideip1) then
            t = fi / (fi - fip1)
            m = m + 1
            if (m > MAXV) stop 'clip_by_point: MAXV exceeded'
            newx(m) = px(i) + t * (px(ip1) - px(i))
            newy(m) = py(i) + t * (py(ip1) - py(i))
         end if
      end do
      n = m
      px(1:n) = newx(1:n)
      py(1:n) = newy(1:n)
   end subroutine clip_by_point

   ! Build one real cell's Voronoi polygon (world coordinates) around seed
   ! point P=(px0,py0), given the full tiled seed-point set (allpx,allpy)
   ! and a spatial grid over it. selfidx is P's own index in allpx/allpy
   ! (excluded as a candidate). Returns polygon vertices in cellx/celly and
   ! their count in ncell.
   subroutine build_cell_polygon(px0, py0, selfidx, allpx, allpy, g, big, &
         cellx, celly, ncell)
      real(dp), intent(in) :: px0, py0, big
      integer, intent(in) :: selfidx
      real(dp), intent(in) :: allpx(:), allpy(:)
      type(grid_t), intent(in) :: g
      real(dp), intent(out) :: cellx(MAXV), celly(MAXV)
      integer, intent(out) :: ncell

      integer, parameter :: MAXCAND = 4096
      integer :: cand(MAXCAND), ncand, i, best, tmpi, nkeep
      integer :: keepidx(MAXCAND)
      real(dp) :: keepdist(MAXCAND), tmpr
      real(dp) :: R, maxR, qx, qy
      logical :: done
      integer :: iter

      ! Gather ALL candidates within a growing radius R, sorted by distance,
      ! clip against every one of them (from scratch each time R grows --
      ! cheap at this problem scale), and accept the result once the
      ! resulting polygon's farthest vertex is provably unreachable by
      ! anything beyond R (maxR <= R/2: see plan/log.txt for the proof).
      R = 4.0d0 * g%cell
      done = .false.
      do iter = 1, 30
         call grid_query(g, px0, py0, R, cand, ncand, MAXCAND)

         nkeep = 0
         do i = 1, ncand
            if (cand(i) == selfidx) cycle
            tmpr = sqrt((allpx(cand(i)) - px0)**2 + (allpy(cand(i)) - py0)**2)
            if (tmpr < 1.0d-12) cycle   ! coincident point safeguard
            if (tmpr > R) cycle         ! grid_query over-returns near bucket edges
            nkeep = nkeep + 1
            keepidx(nkeep) = cand(i)
            keepdist(nkeep) = tmpr
         end do

         ! selection sort by distance ascending (nkeep is small -- a few
         ! dozen to low hundreds typically -- O(nkeep^2) is negligible)
         do i = 1, nkeep - 1
            best = i
            do tmpi = i + 1, nkeep
               if (keepdist(tmpi) < keepdist(best)) best = tmpi
            end do
            if (best /= i) then
               tmpr = keepdist(i); keepdist(i) = keepdist(best); keepdist(best) = tmpr
               tmpi = keepidx(i); keepidx(i) = keepidx(best); keepidx(best) = tmpi
            end if
         end do

         ! initial polygon: a large CCW box, P-centered
         ncell = 4
         cellx(1) = -big; celly(1) = -big
         cellx(2) =  big; celly(2) = -big
         cellx(3) =  big; celly(3) =  big
         cellx(4) = -big; celly(4) =  big

         do i = 1, nkeep
            qx = allpx(keepidx(i)) - px0
            qy = allpy(keepidx(i)) - py0
            call clip_by_point(cellx, celly, ncell, qx, qy)
         end do

         call poly_maxdist(cellx, celly, ncell, maxR)
         if (maxR <= 0.5d0 * R) then
            done = .true.
            exit
         end if
         R = R * 2.0d0
      end do

      if (.not. done) then
         write(*,*) 'build_cell_polygon: candidate radius exhausted for a cell -- increase MAXCAND/search radius'
         stop 1
      end if
   end subroutine build_cell_polygon

   subroutine poly_maxdist(px, py, n, maxR)
      real(dp), intent(in) :: px(MAXV), py(MAXV)
      integer, intent(in) :: n
      real(dp), intent(out) :: maxR
      integer :: i
      real(dp) :: d
      maxR = 0.0d0
      do i = 1, n
         d = sqrt(px(i)**2 + py(i)**2)
         if (d > maxR) maxR = d
      end do
   end subroutine poly_maxdist

   ! Rotate a cell's global-vertex-index list to start at its min-x vertex,
   ! then force clockwise winding via a shoelace-sign check -- identical
   ! convention to StoreData.m's polyorient-based flip. (StoreData.m also
   ! does an intermediate "flip if a(2)>a(end)" pass before this, but that
   ! only affects which of the two possible orientations a rotation starts
   ! from -- something this final clockwise-enforcement step already fully
   ! determines regardless of starting direction, so it is redundant here
   ! and omitted.)
   subroutine reorder_cell(gidx, m, gx, gy)
      integer, intent(inout) :: gidx(MAXV)
      integer, intent(in) :: m
      real(dp), intent(in) :: gx(:), gy(:)
      integer :: i, ip1, kmin, tmp(MAXV)
      real(dp) :: xmin_val, area

      kmin = 1
      xmin_val = gx(gidx(1))
      do i = 2, m
         if (gx(gidx(i)) < xmin_val) then
            xmin_val = gx(gidx(i))
            kmin = i
         end if
      end do
      if (kmin > 1) then
         tmp(1:m) = cshift(gidx(1:m), kmin - 1)
         gidx(1:m) = tmp(1:m)
      end if

      area = 0.0d0
      do i = 1, m
         ip1 = mod(i, m) + 1
         area = area + gx(gidx(i)) * gy(gidx(ip1)) - gy(gidx(i)) * gx(gidx(ip1))
      end do
      area = 0.5d0 * area

      if (area > 0.0d0) then
         tmp(2:m) = gidx(m:2:-1)
         gidx(2:m) = tmp(2:m)
      end if
   end subroutine reorder_cell

end module mesh_gen_utils


program generate_initial_mesh
   use mesh_gen_utils
   use spatial_grid
   implicit none

   integer :: Lx, Ly, seed, cell_headroom, vertex_slot_headroom, vertex_pool_headroom
   real(dp) :: jitter_interior_x, jitter_interior_y, jitter_boundary_x, jitter_boundary_y

   integer :: N, i, j, k, iy_lat, ix_lat
   real(dp), allocatable :: x0(:), y0(:)
   logical, allocatable :: is_boundary(:)
   real(dp) :: offset, r

   real(dp), allocatable :: allpx(:), allpy(:)
   integer :: Ntiled
   type(grid_t) :: seedgrid, vertgrid
   real(dp) :: filter_seed_cell, big

   real(dp), allocatable :: gx(:), gy(:)
   integer :: nglobal, maxglobal
   integer, allocatable :: inn_raw(:,:)
   integer, allocatable :: num_raw(:)
   integer :: maxvert_seen

   real(dp) :: cellx(MAXV), celly(MAXV)
   integer :: gidx(MAXV), ncell_local

   integer, allocatable :: inn(:,:), num(:)
   real(dp), allocatable :: v(:,:)
   integer :: v_dim1, v_dim2, inn_dim1, inn_dim2, num_dim

   integer, parameter :: MAXMATCH = 64
   integer :: gidx_found
   real(dp) :: eps_merge, mergecell

   integer :: seedsz
   integer, allocatable :: seedarr(:)

   ! ---------------- read config ----------------
   open(unit=10, file='para_MeshGen.dat', status='old')
   read(10,*) Lx
   read(10,*) Ly
   read(10,*) seed
   read(10,*) jitter_interior_x
   read(10,*) jitter_interior_y
   read(10,*) jitter_boundary_x
   read(10,*) jitter_boundary_y
   read(10,*) cell_headroom
   read(10,*) vertex_slot_headroom
   read(10,*) vertex_pool_headroom
   close(10)

   if (mod(Ly, 2) /= 0) then
      write(*,*) 'Generate_Initial_Mesh: Ly must be even (hexagonal-offset lattice assumption).'
      stop 1
   end if
   if (Lx < 3 .or. Ly < 3) then
      write(*,*) 'Generate_Initial_Mesh: Lx and Ly must both be >= 3.'
      stop 1
   end if

   N = Lx * Ly

   ! ---------------- seed the RNG ----------------
   call random_seed(size=seedsz)
   allocate(seedarr(seedsz))
   do i = 1, seedsz
      seedarr(i) = seed + 104729 * (i - 1)   ! spread a single int across the seed array
   end do
   call random_seed(put=seedarr)

   ! ---------------- hexagonal-offset lattice + anisotropic jitter ----------------
   allocate(x0(N), y0(N), is_boundary(N))
   do k = 1, N
      ix_lat = (k - 1) / Ly + 1        ! x-column index, 1..Lx
      iy_lat = k - (ix_lat - 1) * Ly   ! y-row index, 1..Ly
      offset = 0.0d0
      if (mod(iy_lat, 2) == 0) offset = 0.5d0
      x0(k) = dble(ix_lat) - offset
      y0(k) = dble(iy_lat)
      is_boundary(k) = (iy_lat == 1 .or. iy_lat == Ly .or. ix_lat == 1 .or. ix_lat == Lx)
   end do

   do k = 1, N
      call random_number(r)
      ! Same random draw scales both x and y offsets, exactly mirroring
      ! Main.m's MeshGenerator (which reuses one `ran` array for both) --
      ! this makes jitter direction correlated (not independent x/y noise).
      if (is_boundary(k)) then
         x0(k) = x0(k) + jitter_boundary_x * (2.0d0 * r - 1.0d0)
         y0(k) = y0(k) + jitter_boundary_y * (2.0d0 * r - 1.0d0)
      else
         x0(k) = x0(k) + jitter_interior_x * (2.0d0 * r - 1.0d0)
         y0(k) = y0(k) + jitter_interior_y * (2.0d0 * r - 1.0d0)
      end if
   end do

   ! ---------------- 8 periodic ghost copies ----------------
   Ntiled = 9 * N
   allocate(allpx(Ntiled), allpy(Ntiled))
   allpx(1:N) = x0; allpy(1:N) = y0
   allpx(N+1:2*N)   = x0 + Lx;      allpy(N+1:2*N)   = y0
   allpx(2*N+1:3*N) = x0 - Lx;      allpy(2*N+1:3*N) = y0
   allpx(3*N+1:4*N) = x0;           allpy(3*N+1:4*N) = y0 + Ly
   allpx(4*N+1:5*N) = x0;           allpy(4*N+1:5*N) = y0 - Ly
   allpx(5*N+1:6*N) = x0 + Lx;      allpy(5*N+1:6*N) = y0 + Ly
   allpx(6*N+1:7*N) = x0 - Lx;      allpy(6*N+1:7*N) = y0 - Ly
   allpx(7*N+1:8*N) = x0 + Lx;      allpy(7*N+1:8*N) = y0 - Ly
   allpx(8*N+1:9*N) = x0 - Lx;      allpy(8*N+1:9*N) = y0 + Ly

   ! spatial grid over the tiled seed points, for candidate gathering.
   ! Bucket size set from jitter magnitudes: true Voronoi neighbors sit
   ! within roughly one lattice spacing plus jitter, so a bucket a few
   ! times the max jitter is generous without being coarse.
   filter_seed_cell = max(0.5d0, 2.0d0 * max(jitter_interior_x, jitter_interior_y, &
        jitter_boundary_x, jitter_boundary_y))
   call grid_init(seedgrid, -2.0d0*Lx, -2.0d0*Ly, 3.0d0*Lx, 3.0d0*Ly, filter_seed_cell, Ntiled)
   do i = 1, Ntiled
      call grid_insert(seedgrid, allpx(i), allpy(i), i)
   end do

   big = 10.0d0 * dble(max(Lx, Ly))

   ! ---------------- build each real cell's polygon ----------------
   maxglobal = 8 * N + 512
   allocate(gx(maxglobal), gy(maxglobal))
   nglobal = 0

   ! merge grid over the (growing) global vertex list
   eps_merge = 1.0d-6
   mergecell = max(50.0d0 * eps_merge, 0.02d0)
   call grid_init(vertgrid, -1.0d0, -1.0d0, dble(Lx) + 1.0d0, dble(Ly) + 1.0d0, mergecell, maxglobal)

   allocate(inn_raw(N, MAXV), num_raw(N))
   inn_raw = 0
   maxvert_seen = 0

   do k = 1, N
      call build_cell_polygon(x0(k), y0(k), k, allpx, allpy, seedgrid, big, cellx, celly, ncell_local)

      do i = 1, ncell_local
         ! shift back to world coordinates before merging
         call find_or_add_vertex(vertgrid, gx, gy, nglobal, maxglobal, &
              cellx(i) + x0(k), celly(i) + y0(k), eps_merge, gidx_found)
         gidx(i) = gidx_found
      end do

      call reorder_cell(gidx, ncell_local, gx, gy)

      inn_raw(k, 1:ncell_local) = gidx(1:ncell_local)
      num_raw(k) = ncell_local
      if (ncell_local > maxvert_seen) maxvert_seen = ncell_local
   end do

   ! ---------------- assemble final padded arrays ----------------
   num_dim = N + cell_headroom
   inn_dim2 = num_dim
   inn_dim1 = maxvert_seen + vertex_slot_headroom
   v_dim1 = 2
   v_dim2 = nglobal + vertex_pool_headroom

   allocate(inn(inn_dim1, inn_dim2)); inn = 0
   allocate(num(num_dim)); num = 0
   allocate(v(v_dim1, v_dim2)); v = 0.0d0

   do k = 1, N
      num(k) = num_raw(k)
      inn(1:num_raw(k), k) = inn_raw(k, 1:num_raw(k))
   end do
   do j = 1, nglobal
      v(1, j) = gx(j)
      v(2, j) = gy(j)
   end do

   ! ---------------- write output files ----------------
   open(unit=20, file='v_in.dat', status='replace')
   do j = 1, v_dim2
      write(20,*) v(1, j), v(2, j)
   end do
   close(20)

   open(unit=21, file='inn_in.dat', status='replace')
   do k = 1, inn_dim2
      write(21,*) (inn(j, k), j = 1, inn_dim1)
   end do
   close(21)

   open(unit=22, file='num_in.dat', status='replace')
   do k = 1, num_dim
      write(22,*) num(k)
   end do
   close(22)

   open(unit=23, file='para_MeshDims.dat', status='replace')
   write(23,*) Lx
   write(23,*) Ly
   write(23,*) num_dim
   write(23,*) v_dim1
   write(23,*) v_dim2
   write(23,*) inn_dim1
   write(23,*) inn_dim2
   close(23)

   write(*,*) 'Generate_Initial_Mesh: wrote', N, 'cells,', nglobal, 'vertices.'
   write(*,*) 'para_MeshDims.dat: Lx,Ly,num_dim,v_dim1,v_dim2,inn_dim1,inn_dim2 =', &
        Lx, Ly, num_dim, v_dim1, v_dim2, inn_dim1, inn_dim2

contains

   ! Look up a world-coordinate vertex in the global list within eps; add it
   ! if not already present. Grows gx/gy in place (both sized maxn) and
   ! inserts into the shared merge grid.
   subroutine find_or_add_vertex(g, gx, gy, nglobal, maxn, x, y, eps, idx_out)
      type(grid_t), intent(inout) :: g
      real(dp), intent(inout) :: gx(:), gy(:)
      integer, intent(inout) :: nglobal
      integer, intent(in) :: maxn
      real(dp), intent(in) :: x, y, eps
      integer, intent(out) :: idx_out
      integer :: cand(MAXMATCH), ncand, i
      real(dp) :: d2, eps2

      eps2 = eps * eps
      call grid_query(g, x, y, g%cell, cand, ncand, MAXMATCH)
      do i = 1, ncand
         d2 = (gx(cand(i)) - x)**2 + (gy(cand(i)) - y)**2
         if (d2 <= eps2) then
            idx_out = cand(i)
            return
         end if
      end do

      nglobal = nglobal + 1
      if (nglobal > maxn) stop 'find_or_add_vertex: increase maxglobal (vertex_pool headroom)'
      gx(nglobal) = x
      gy(nglobal) = y
      call grid_insert(g, x, y, nglobal)
      idx_out = nglobal
   end subroutine find_or_add_vertex

end program generate_initial_mesh
