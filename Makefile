ifeq ($(origin FC),default)
FC = mpifort
endif
FFLAGS  ?= -O3 -march=native
NP      ?= 4
MPIRUN  ?= mpirun
GMSH    ?= gmsh

ifeq ($(DEBUG),1)
FFLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall
endif
FFLAGS += -std=f2018 -Jbuild

SRC = kinds sort mesh partition config cases domain solver output main
OBJ = $(SRC:%=build/%.o)

MESHES = meshes/bump.msh meshes/thacker.msh meshes/island.msh meshes/island_fine.msh \
         meshes/dambreak_1.msh meshes/dambreak_2.msh meshes/dambreak_3.msh meshes/dambreak_4.msh

.PHONY: all meshes test clean

all: estran

estran: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

build/%.o: src/%.f90 | build
	$(FC) $(FFLAGS) -c $< -o $@

build:
	mkdir -p build

# modules must be compiled before their users
build/mesh.o: build/kinds.o build/sort.o
build/partition.o: build/kinds.o build/sort.o
build/config.o: build/kinds.o
build/cases.o: build/config.o
build/domain.o: build/mesh.o build/partition.o
build/solver.o: build/domain.o build/cases.o
build/output.o: build/solver.o
build/main.o: build/output.o

meshes: $(MESHES)

meshes/%.msh: geo/%.geo
	@mkdir -p meshes
	$(GMSH) -2 -format msh22 -v 2 $< -o $@

# the dam-break sequence halves the mesh size each time
meshes/dambreak_%.msh: geo/dambreak.geo
	@mkdir -p meshes
	$(GMSH) -2 -format msh22 -v 2 -setnumber lc $$(python3 -c "print(1.0/2**$*)") $< -o $@

meshes/island_fine.msh: geo/island.geo
	@mkdir -p meshes
	$(GMSH) -2 -format msh22 -v 2 -setnumber li 0.035 $< -o $@

test: estran meshes
	MPIRUN="$(MPIRUN)" NP=$(NP) ./tests/run_tests.sh

clean:
	rm -rf build estran
