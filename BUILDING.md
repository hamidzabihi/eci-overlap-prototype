# Building wfoverlap with gfortran

The SHARC `wfoverlap` Makefile defaults to Intel `ifx`. This documents
how to build it with `gfortran` on a standard Linux system.

## Prerequisites

    sudo apt install gfortran liblapack-dev libblas-dev

## Steps

    cd ~/sharc4/wfoverlap/source

    # Back up the Intel-specific readers, use the ASCII/dummy variants
    mv read_dalton.f90 read_dalton.f90.bak
    mv read_molcas.f90 read_molcas.f90.bak
    cp read_dalton_dummy.f90 read_dalton.f90
    cp read_molcas_dummy.f90 read_molcas.f90

    # Clean and build
    rm -f *.o *.mod *.x *_genmod.f90
    make FC=gfortran \
         FCFLAGS="-O3 -fopenmp -cpp -fdefault-integer-8 -DEXTBLAS -ffree-line-length-none -fbacktrace" \
         LINKFLAGS="-fopenmp" \
         LIBS="-llapack -lblas"

    # Install
    cp wfoverlap.x ../../bin/

## Notes

- `-fdefault-integer-8` is essential; without it, the determinant
  routines silently produce wrong results.
- The dummy readers (`read_dalton_dummy.f90`, `read_molcas_dummy.f90`)
  are needed for a standalone build; they provide stub implementations
  of the Molcas/Dalton interfaces that are only invoked by the native
  ONEINT/aoints formats, not by the ASCII path.
