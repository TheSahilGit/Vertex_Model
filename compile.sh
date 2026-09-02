gfortran -pg -O3 -march=native -fbounds-check -fbacktrace allocation.f90 array_info.f90 Geometry.f90  T1_swap.f90 T2_swap.f90  System_Info.f90 Force.f90  Stress.f90 Proliferation.f90 vertexmain.f90 -o  vertexmain.exe

# Standalone initial-mesh generator (replaces Main.m -- see log.txt); reads
# para_MeshGen.dat, writes v_in.dat/inn_in.dat/num_in.dat/para_MeshDims.dat.
# Independent of the files above, so it's a separate build line.
# -->
# Only compile for nrun = 1, otherwise rewrites para_MeshDims.dat
# -->
#gfortran -O3 -fbounds-check -fbacktrace Generate_Initial_Mesh.f90 -o generate_initial_mesh.exe
