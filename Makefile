# Vertex Model build system.
#
#   make                 -> builds bin/vertexmain, WITHOUT -fbounds-check/
#                            -fbacktrace (fast -- use for production runs)
#   make fbounds         -> builds bin/vertexmain WITH -fbounds-check/
#                            -fbacktrace (slower -- array-bounds checking
#                            has real runtime cost; use while debugging,
#                            not for production runs)
#   make initial_config  -> builds bin/generate_initial_mesh
#                            (standalone mesh generator -- only needed to
#                            (re)create mesh/v_in.dat/mesh/inn_in.dat/
#                            mesh/num_in.dat/mesh/para_MeshDims.dat; see
#                            PARAMETERS.md)
#   make clean           -> removes build/ and bin/
#   make dataclean       -> removes simulation output from data/, keeping
#                            only the small tracked reference snapshot
#                            (the four data/*00050000*/motility_store.dat
#                            files .gitignore already carves out) -- see
#                            below for exactly what this does and why.
#
# Both binaries must be RUN FROM THE REPO ROOT, not from bin/ -- they open
# their input/output files (para_Simulation.dat, para_MeshGen.dat,
# mesh/..., data/...) via relative paths:
#   ./bin/vertexmain
#   ./bin/generate_initial_mesh
#
# Compiler/flags/source order are unchanged from the old compile.sh (see
# legacy/compile.sh) -- only file locations, the mesh/ input paths, and
# (this entry) the optional -fbounds-check/-fbacktrace toggle have changed
# since this Makefile was first added (log.txt).

FC = gfortran
FFLAGS_BASE  = -pg -O3 -march=native
FFLAGS_DEBUG = -fbounds-check -fbacktrace
FFLAGS_MESH  = -O3 -fbounds-check -fbacktrace

FFLAGS_MAIN_default = $(FFLAGS_BASE)
FFLAGS_MAIN_fbounds  = $(FFLAGS_BASE) $(FFLAGS_DEBUG)

SRC_DIR   = src
BUILD_DIR = build
BIN_DIR   = bin
DATA_DIR  = data

# Dependency order matters (gfortran compiles a multi-file command in the
# order given; a module must appear before any file that `use`s it) --
# identical order to the original compile.sh.
VERTEXMAIN_SRCS = $(SRC_DIR)/allocation.f90 \
                   $(SRC_DIR)/array_info.f90 \
                   $(SRC_DIR)/Geometry.f90 \
                   $(SRC_DIR)/T1_swap.f90 \
                   $(SRC_DIR)/T2_swap.f90 \
                   $(SRC_DIR)/System_Info.f90 \
                   $(SRC_DIR)/Force.f90 \
                   $(SRC_DIR)/Stress.f90 \
                   $(SRC_DIR)/Proliferation.f90 \
                   $(SRC_DIR)/vertexmain.f90

.PHONY: all fbounds initial_config clean dataclean

all: $(BIN_DIR)/vertexmain

fbounds: $(VERTEXMAIN_SRCS) $(BUILD_DIR)/.flags-fbounds | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(FFLAGS_MAIN_fbounds) -J$(BUILD_DIR) $(VERTEXMAIN_SRCS) -o $(BIN_DIR)/vertexmain

initial_config: $(BIN_DIR)/generate_initial_mesh

# A per-variant sentinel file, so switching between `make` and
# `make fbounds` forces vertexmain to relink with the right flags even
# though no .f90 source changed -- without it, a binary built by one
# variant would look "up to date" to make when the other variant is
# requested next, and silently NOT be rebuilt (log.txt: target-specific
# variables were tried first and don't propagate into a prerequisite
# NAME the way this needs, hence the two explicit recipes below rather
# than one parameterized rule).
$(BIN_DIR)/vertexmain: $(VERTEXMAIN_SRCS) $(BUILD_DIR)/.flags-default | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(FFLAGS_MAIN_default) -J$(BUILD_DIR) $(VERTEXMAIN_SRCS) -o $@

$(BUILD_DIR)/.flags-default: | $(BIN_DIR) $(BUILD_DIR)
	@rm -f $(BUILD_DIR)/.flags-fbounds
	@touch $@

$(BUILD_DIR)/.flags-fbounds: | $(BIN_DIR) $(BUILD_DIR)
	@rm -f $(BUILD_DIR)/.flags-default
	@touch $@

$(BIN_DIR)/generate_initial_mesh: $(SRC_DIR)/Generate_Initial_Mesh.f90 | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(FFLAGS_MESH) -J$(BUILD_DIR) $< -o $@

$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f gmon.out

# Uses git itself to know what's regenerable vs. kept, rather than
# hardcoding filenames here a second time: `git clean -X` removes only
# files git already considers ignored (per .gitignore's `data/*` +
# `!data/...` carve-outs), so the four tracked reference files always
# survive, and this stays correct automatically if those carve-outs ever
# change. Requires data/ to actually be inside this git repo.
dataclean:
	git clean -X -f -d -- $(DATA_DIR)
