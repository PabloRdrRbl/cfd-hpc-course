# This is a Makefile, recipes are annoted here and made with make

.phony: all, clean, tidy

PROGRAM_NAMES= test_interpolation potential_flow_solver
PROGRAMS     = $(addprefix bin/,$(PROGRAM_NAMES))

all:    $(PROGRAMS)

clean:
	rm -f tmp/*.mod tmp/*.o

tidy:   clean
	rm -f $(PROGRAMS)

# Compiler & flags -------------------------------------------------------------

FC     := mpifort
FFLAGS := -O3 -ffree-form -fopenmp -Wall -std=gnu -J tmp/ -I tmp/

# Rules ------------------------------------------------------------------------

# generate binary, use all available object files
bin/%: tmp/%.o tmp/kind_parameters.o tmp/constants.o
	$(FC) $(FFLAGS)   $^ -o $@

# generate object files from source files
tmp/%.o: src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Module dependencies ----------------------------------------------------------

tmp/constants.o:                                                               \
        tmp/kind_parameters.o

tmp/spectral_element_interpolation.o:                                          \
        tmp/kind_parameters.o
tmp/test_interpolation.o:                                                      \
        tmp/kind_parameters.o                                                  \
        tmp/spectral_element_interpolation.o

tmp/gauss_jacobi.o:                                                            \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o

tmp/helmholtz_2d.o:                                                            \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o                                                        \
        tmp/gauss_jacobi.o

tmp/potential_flow_solver.o:                                                   \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o                                                        \
        tmp/helmholtz_2d.o

# Program dependencies .........................................................

bin/test_interpolation:                                                        \
        tmp/spectral_element_interpolation.o

bin/potential_flow_solver:                                                     \
	tmp/gauss_jacobi.o                                                     \
        tmp/helmholtz_2d.o
