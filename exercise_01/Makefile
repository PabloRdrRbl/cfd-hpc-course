# This is a Makefile, recipes are annoted here and made with make

.phony: all, clean, tidy

PROGRAM_NAMES= stencil_test                                                    \
               potpourri_a1                                                    \
               potpourri_a2                                                    \
               potpourri_a3                                                    \
               matrix_matrix_test

PROGRAMS     = $(addprefix bin/,$(PROGRAM_NAMES))

all:    $(PROGRAMS)

clean:
	rm -f *.mod *.o tmp/*.mod tmp/*.o

tidy:   clean
	rm -f $(PROGRAMS)

# Compiler & flags -------------------------------------------------------------

FC     := mpifort
FFLAGS := -O3 -ffree-form -Wall -Wpedantic -J tmp/ -I tmp/

# Rules ------------------------------------------------------------------------

# generate binary, use all available object files
bin/%: tmp/%.o tmp/kind_parameters.o tmp/constants.o
	$(FC) $(FFLAGS)   $^ -o $@

# generate object files from source files
tmp/%.o: src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies -----------------------------------------------------------------

# modules
tmp/constants.o:             tmp/kind_parameters.o

tmp/potpourri_print_mlups.o: tmp/kind_parameters.o

tmp/stencil_test.o:          tmp/kind_parameters.o                             \
                             tmp/constants.o                                   \
                             tmp/stencil_baseline.o                            \
                             tmp/stencil_direct.o                              \
                             tmp/stencil_blocked.o

tmp/stencil_baseline.o:      tmp/kind_parameters.o
tmp/stencil_direct.o:        tmp/kind_parameters.o
tmp/stencil_blocked.o:       tmp/kind_parameters.o

tmp/potpourri_print_mlups.o: tmp/kind_parameters.o
tmp/potpourri_a1.o:          tmp/kind_parameters.o                             \
                             tmp/potpourri_print_mlups.o
tmp/potpourri_a2_derivatives.o:                                                \
                             tmp/kind_parameters.o
tmp/potpourri_a2.o:          tmp/kind_parameters.o                             \
                             tmp/potpourri_print_mlups.o                       \
                             tmp/potpourri_a2_derivatives.o
tmp/potpourri_a3.o:          tmp/kind_parameters.o                             \
                             tmp/potpourri_print_mlups.o

tmp/matrix_matrix_simple.o:  tmp/kind_parameters.o
tmp/matrix_matrix_blocked.o: tmp/kind_parameters.o
tmp/matrix_matrix_test.o:    tmp/kind_parameters.o                             \
                             tmp/constants.o                                   \
                             tmp/matrix_matrix_simple.o                        \
                             tmp/matrix_matrix_blocked.o

# programs

bin/stencil_test:            tmp/kind_parameters.o tmp/constants.o             \
                             tmp/stencil_baseline.o                            \
                             tmp/stencil_direct.o                              \
                             tmp/stencil_blocked.o
bin/potpourri_a1:            tmp/potpourri_print_mlups.o
bin/potpourri_a2:            tmp/potpourri_print_mlups.o                       \
                             tmp/potpourri_a2_derivatives.o
bin/potpourri_a3:            tmp/potpourri_print_mlups.o

bin/matrix_matrix_test:      tmp/matrix_matrix_simple.o                        \
                             tmp/matrix_matrix_blocked.o

# Extra flags ------------------------------------------------------------------
tmp/stencil_baseline.o: FFLAGS +=
tmp/potpourri_%.o:      FFLAGS += -O3
