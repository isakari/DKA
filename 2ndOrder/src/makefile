TARGET = ../bin/a.out

FC = gfortran
USE_OPENMP = 1
# OPT = 1
DBG = 1
OBJS = module_globals.o \
	module_autodiff.o \
	module_dka.o \
	main.o
GFORTRAN = gfortran
IFORT = ifort

ifeq ($(FC), $(GFORTRAN))
FDEBUGS = -g -fbounds-check -Wuninitialized -fbacktrace -Wunused
FFLAGS = -cpp
endif
ifeq ($(FC), $(IFORT))
FDEBUGS = -check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback -g
FFLAGS = -cpp -DIFORT
endif
ifdef DBG
FFLAGS += $(FDEBUGS) -DDBG
endif

CC = gcc
CDEBUGS = -g -fbounds-check -Wuninitialized -Wunused
CFLAGS = -MMD -DUSE_COMPLEX -funroll-loops -funswitch-loops 
ifdef DBG
CFLAGS += $(CDEBUGS)
endif

ifdef OPT
ifeq ($(FC), $(GFORTRAN))
FFLAGS += -O3
LDFLAGS += -O3
CFLAGS += -march=native -O3
endif
ifeq ($(FC), $(IFORT))
FFLAGS += -O3 -ip
LDFLAGS += -O3 -ip
CFLAGS += -march=native -O3
endif
endif

ifdef DBG
LIBS = -lm
endif
ifdef OPT
LIBS = -lm
endif
LIBS += -lslatec
ifeq ($(FC), $(GFORTRAN))
LIBS += -llapack -lblas
endif
ifeq ($(FC), $(IFORT))
LIBS += -mkl
endif

LDFLAGS = $(LIBS)
INCLUDES = 

ifdef USE_OPENMP
ifeq ($(FC), $(GFORTRAN))
CFLAGS += -DUSE_OPENMP -fopenmp
FFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
ifeq ($(FC), $(IFORT))
CFLAGS += -DUSE_OPENMP -fopenmp
FFLAGS += -qopenmp
LDFLAGS += -qopenmp
endif
endif

.SUFFIXES:
.SUFFIXES: .o .f90 .f .c

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(FC) $(OBJS) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod *.d *~ $(TARGET) *genmod*

# dependency
-include *.d
