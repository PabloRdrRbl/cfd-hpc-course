# -— Find files —-

# Find all source files, create a list of corresponding object files
SRCS=$(wildcard *.f)
OBJS=$(patsubst %.f,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
#MODS=$(wildcard mod*.f)
#MOD_OBJS=$(patsubst %.f,%.o,$(MODS))


# -- Compiler configuration --

# Compiler
FC = gfortran
# Compile flags
FFLAGS = -ffree-form -Wall
# Linker flags
LFLAGS = -ffree-form -Wall


# -- Varibles --

PROGRAM_NAMES = hello calculate_e calculate_pi jacobi_iteration
PROGRAMS = $(addprefix bin/,$(PROGRAM_NAMES))

PROGRAM = hello
PRG_OBJ = $(PROGRAM).o

# -- Main targets --

# Use .PHONY to mark targets that don’t correspond to files
.PHONY: default

# make without parameters will make first target found
default : $(PROGRAM)


# -- Compilation rules --

# Clean the suffixes
.SUFFIXES:

# Set the suffixes we are interested in
.SUFFIXES: .f .o

# Compiler steps for all objects
$(OBJS) : %.o : %.f
    $(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
    $(FC) $(FLFLAGS) -o $@ $^

# -- Program dependencies --

tmp/constants.o:	tmp/kind_parameters.o
	
tmp/hello.o:	tmp/constants.o

tmp/jacobi_iteration.o:	tmp/kind_parameters.o tmp/constants.o



# -- Cleaning --

.PHONY:	clean, tidy

clean:
	rm -f *.{mod,o} work/*.{mod,o} rm -rf tmp/*.{mod,o}

tidy:	clean
	rm -f $(PROGRAMS)


# -- Some Make tips --

# Write rules as 'target: dependencies'
# Specify update actions in a tab-indented block under the rule

# Use '$@' to refer to the target of the current rule
# Use '$^' to refer to the dependencies of the current rule
# Use '$<' to refer to the first dependency of the current rule

# Use the wildcard '%' as a placeholder in targets and dependencies
# Use the special variable '$*' to refer to matching sets of files in actions
# Reference variables using '$(...)'

# Use 'wildcard' function to get lists of files matching a pattern
# Use 'patsubst'function to rewrite file names
