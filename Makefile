#Makefile for MCMC code

# gfortran is default compiler
COMPILER ?= gnu

# Compiler-specific settings
ifeq ($(COMPILER),gnu)
 FC ?= gfortran
 CC ?= gcc
 LD := $(FC)

 # Flag for setting default real/double kind.
 REAL8_FLAGS := -fdefault-real-8 -fdefault-double-8

 # Use OPENBLAS?
 OPENBLAS ?= TRUE


ifeq ($(DEBUG),TRUE)
 # Flags that help to diagnose issues
 FCFLAGS := -O0 -g -fbacktrace -Wall -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow
 ifneq ($(OLD_GCC), TRUE)
  FCFLAGS := $(FCFLAGS) -fcheck=no-array-temps
 endif
 # Don't use warnings, pedantic or otherwise, for bin model code.
 FCFLAGS_BIN := -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
 CFLAGS_RANDNUM := -std=gnu99
else
 # Flags for performance
 FCFLAGS := -O2 -march=native
 FLFLAGS_BIN := $(FCFLAGS)
 # RandNum settings shamelessly copied from benchmark Makefile.
 FCFLAGS_RANDNUM := -Ofast -march=native
 CFLAGS_RANDNUM := -Ofast -march=native -std=gnu99
 CPPFLAGS_RANDNUM := -DHAVE_SSE2
endif # DEBUG

ifeq ($(USE_OPENMP),TRUE)
  FCFLAGS := $(FCFLAGS) -fopenmp
  LDFLAGS := $(LDFLAGS) -fopenmp
endif
$(info    FCFLAGS is $(FCFLAGS))

else ifeq ($(COMPILER),intel)
 FC ?= ifort
 CC ?= icc
 LD := $(FC)
 LDFLAGS += -mkl

 # Flag for setting default real kind.
 REAL8_FLAGS = -r8

 # Use OPENBLAS?
 OPENBLAS ?= FALSE

ifeq ($(DEBUG),TRUE)
 # Flags that help to diagnose issues
 FCFLAGS := -O0 -g -traceback -stand f18 -warn all -check all,noarg_temp_created -fpe0
 # Turn off warnings for bin code.
 FCFLAGS_BIN := -O0 -g -traceback -check all,noarg_temp_created -fpe0

 CFLAGS_RANDNUM := -std=c99
else
 # Flags for performance
 FCFLAGS := -O2 -xHost
 FCFLAGS_BIN := $(FCFLAGS)
 # RandNum settings shamelessly copied from benchmark Makefile.
 FFLAGS_RANDNUM := -O3 -xHost -fp-model fast -mkl -no-prec-div -no-prec-sqrt -override-limits
 CFLAGS_RANDNUM := -O3 -xHost -fp-model fast -std=c99
 CPPDEFS_RANDNUM = -DINTEL_MKL -DHAVE_SSE2
endif # DEBUG

ifeq ($(USE_OPENMP),TRUE)
 FCFLAGS := $(FCFLAGS) -qopenmp
 LDFLAGS := $(LDFLAGS) -qopenmp
endif # USEOPENMP
#$(info $$FCFLAGS is [${FCFLAGS}])
$(info    FCFLAGS is $(FCFLAGS))


endif # COMPILER

ifeq ($(POSTERIOR),TRUE)
  MAINFILE := posterior_main
else
  MAINFILE := main
endif 

# Power used for double precision SIMD Fast Mersenne Twister.
CPPFLAGS_RANDNUM := $(CPPFLAGS_RANDNUM) -DDSFMT_MEXP=19937

# Autopromote all real/double intrinsics to 8 bytes.
FCFLAGS := $(FCFLAGS) $(REAL8_FLAGS)
FCFLAGS_BIN := $(FCFLAGS_BIN) $(REAL8_FLAGS)

# Location of the netCDF installation.
NETCDF ?= /usr/local

# CIME shared code directory
CIMESHR_DIR := ./cime_share

vpath %.F90 $(CIMESHR_DIR)

# Random number generator code directory
RANDNUM_DIR := $(CIMESHR_DIR)/RandNum

# Include flags
INCS := -I$(NETCDF)/include -I./$(RANDNUM_DIR)/include

# External libraries for linking
LDLIBS := -lnetcdff -lnetcdf
ifeq ($(OPENBLAS),TRUE)
 LDLIBS := $(LDLIBS) -lopenblas
endif

# Linking flags
LDFLAGS := $(LDFLAGS) -L$(NETCDF)/lib

.DEFAULT:
	-echo $@ does not exist.
.PHONY: clean all
.SUFFIXES: .F90 .f90 .c .o

all: mcmc.x

.F90.o:
	$(FC) -c $(INCS) $(CPPFLAGS) $(FCFLAGS) $<

.f90.o:
	$(FC) -c $(INCS) $(FCFLAGS) $<

.c.o:
	$(CC) -c $(INCS) $(CPPFLAGS) $(CFLAGS) $<

%.o: %.mod

# RandNum file compilation.
dSFMT.o: $(RANDNUM_DIR)/src/dsfmt_f03/dSFMT.c
	$(CC) $(INCS) $(CPPFLAGS_RANDNUM) $(CFLAGS_RANDNUM) -c $<
dSFMT_interface.o: $(RANDNUM_DIR)/src/dsfmt_f03/dSFMT_interface.F90
	$(FC) $(INCS) $(CPPFLAGS_RANDNUM) $(FCFLAGS_RANDNUM) -c $<
dSFMT_utils.o: $(RANDNUM_DIR)/src/dsfmt_f03/dSFMT_utils.c $(RANDNUM_DIR)/include/dSFMT.h
	$(CC) $(INCS) $(CPPFLAGS_RANDNUM) $(CFLAGS_RANDNUM) -c $<
kissvec_mod.o: $(RANDNUM_DIR)/src/kissvec/kissvec_mod.F90
	$(FC) $(INCS) $(CPPFLAGS_RANDNUM) $(FCFLAGS_RANDNUM) -c $<
kissvec.o: $(RANDNUM_DIR)/src/kissvec/kissvec.c
	$(CC) $(INCS) $(CPPFLAGS_RANDNUM) $(CFLAGS_RANDNUM) -c $<
mersennetwister_mod.o: $(RANDNUM_DIR)/src/mt19937/mersennetwister_mod.F90
	$(FC) $(INCS) $(CPPFLAGS_RANDNUM) $(FCFLAGS_RANDNUM) -c $<
shr_RandNum_mod.o: $(RANDNUM_DIR)/src/shr_RandNum_mod.F90 kissvec_mod.o mersennetwister_mod.o dSFMT_interface.o
	$(FC) $(INCS) $(CPPFLAGS_RANDNUM) $(FCFLAGS_RANDNUM) -c $<

OBJS_SHARE := shr_kind_mod.o shr_const_mod.o shr_wv_sat_mod.o
OBJS_RANDNUM := kissvec_mod.o mersennetwister_mod.o dSFMT.o dSFMT_interface.o \
                dSFMT_utils.o shr_RandNum_mod.o kissvec.o
OBJS_UTIL := io_util.o netcdf_util.o rand_util.o micro_util.o
OBJS_MICRO := micro_driver.o simple_boss.o morr_bulk.o mp_bin.o module_mp_bin_habit.o advect1d.o

OBJS_OFFLINE := offline_micro.o $(OBJS_MICRO) $(OBJS_SHARE) $(OBJS_UTIL) $(OBJS_RANDNUM)

OBJS_MCMC := $(MAINFILE).o fwd_model.o rand_util.o $(OBJS_MICRO) $(OBJS_SHARE) \
	$(OBJS_RANDNUM) $(OBJS_UTIL)

# Bin model compilation
module_mp_bin_habit.o: module_mp_bin_habit.f90 $(OBJS_SHARE)
	$(FC) -c $(INCS) $(FCFLAGS_BIN) $<

# Dependency of Fortran modules on one another.
shr_const_mod.o: shr_kind_mod.o
rand_util.o: $(OBJS_SHARE) $(OBJS_RANDNUM)
netcdf_util.o: $(OBJS_SHARE) io_util.o
micro_util.o: $(OBJS_SHARE)
advect1d.o: $(OBJS_SHARE)
morr_bulk.o: $(OBJS_SHARE) $(OBJS_UTIL)
mp_bin.o: $(OBJS_SHARE) $(OBJS_UTIL) module_mp_bin_habit.o
simple_boss.o: $(OBJS_SHARE) $(OBJS_UTIL)
micro_driver.o: $(OBJS_SHARE) morr_bulk.o simple_boss.o mp_bin.o advect1d.o $(OBJS_UTIL)
fwd_model.o: $(OBJS_MICRO) $(OBJS_SHARE) $(OBJS_UTIL)
$(MAINFILE).o: fwd_model.o $(OBJS_SHARE) $(OBJS_UTIL)

# Offline driver
offline_micro.o: $(OBJS_MICRO) $(OBJS_SHARE)

# Linking stage
mcmc.x: $(OBJS_MCMC)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Offline driver linking
offline_micro.x: $(OBJS_OFFLINE)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

#..clean
clean :
	rm *.o *.mod
