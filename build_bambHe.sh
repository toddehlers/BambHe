#!/bin/bash
#
# build_direct.sh - compiles forward model executable
#
# Portland Group Fortran compiler
#pgf90 -O3 -c fortran_src_code/*.f fortran_src_code/*.f90
#pgf90 -O3 -o bambHe *.o

# Portland Group Fortran compiler safe
#pgf90 -c fortran_src_code/*.f fortran_src_code/*.f90
#pgf90 -o bambHe *.o

# Portland Group Fortran compiler w/ debug/warnings
#pgf90 -C -w -c fortran_src_code/*.f fortran_src_code/*.f90
#pgf90 -C -w -o bambHe *.o

# Intel Fortran compiler
ifort -O3 -c src/*_module.f
ifort -O3 -c src/*.f src/*.f90
ifort -O3 -o bambHe *.o

# Intel Fortran compiler optimized for Core 2 Duos
#ifort -fast -c fortran_src_code/*.f fortran_src_code/*.f90
#ifort -fast -o bambHe *.o

# Intel Fortran compiler w/ debug/warning flags
#ifort -c -C -w -debug all fortran_src_code/*.f fortran_src_code/*.f90
#ifort -o bambHe -C -w -debug all *.o

# Intel Fortran compiler for idb debugger
#ifort -c -g -pg fortran_src_code/*.f fortran_src_code/*.f90
#ifort -o bambHe -g -pg *.o

# Latest GNU Fortran compiler (slower than ifort or pgf90, but free)
#gfortran -O3 -c fortran_src_code/*.f fortran_src_code/*.f90
#gfortran -O3 -o bambHe *.o

# Latest GNU Fortran compiler w/ debug/warnings
#gfortran -C -w -c fortran_src_code/*.f fortran_src_code/*.f90
#gfortran -C -w -o bambHe *.o

# Another open source Fortran compiler
#g95 -c fortran_src_code/*.f fortran_src_code/*.f90
#g95 -o bambHe *.o

# Remove object files after compilation
rm -f *.o
rm -f *.mod
