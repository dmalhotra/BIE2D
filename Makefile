SCTL_INCLUDE_DIR = SCTL/include

#CXX=c++ # requires g++-8 or newer / icpc (with gcc compatibility 7.5 or newer) / clang++ with llvm-10 or newer
CXX=mpicxx -DSCTL_HAVE_MPI
#CXX=/mnt/sw/nix/store/jrw0k2lr4i16pn5ja1rp34wzazdq7ivw-intel-oneapi-compilers-2023.0.0/compiler/2023.0.0/linux/bin/icpx
CXXFLAGS = -std=c++11 -fopenmp -Wall -Wfloat-conversion # need C++11 and OpenMP

#Optional flags
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS += -O0 -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all -fstack-protector # debug build
  CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
else
	CXXFLAGS += -O3 -march=native -DNDEBUG # release build
endif

OS = $(shell uname -s)
ifeq "$(OS)" "Darwin"
	CXXFLAGS += -g -rdynamic -Wl,-no_pie # for stack trace (on Mac)
else
	CXXFLAGS += -gdwarf-4 -g -rdynamic # for stack trace -gstrict-dwarf
endif

CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling

CXXFLAGS += -DSCTL_QUAD_T=__float128 # Enable quadruple precision

#CXXFLAGS += -lblas -DSCTL_HAVE_BLAS # use BLAS
#CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK # use LAPACK
#CXXFLAGS += -qmkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK -DSCTL_HAVE_FFTW3_MKL # use MKL BLAS, LAPACK and FFTW (Intel compiler)
#CXXFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK (non-Intel compiler)
#CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK
CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK

#CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
#CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
#CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL

#CXXFLAGS += -DSCTL_HAVE_SVML

#CXXFLAGS += -I${PETSC_DIR}/include -I${PETSC_DIR}/../include -DSCTL_HAVE_PETSC
#LDLIBS += -L${PETSC_DIR}/lib -lpetsc

#PVFMM_INC_DIR = ../include
#PVFMM_LIB_DIR = ../lib/.libs
#CXXFLAGS += -DSCTL_HAVE_PVFMM -I$(PVFMM_INC_DIR)
#LDLIBS += $(PVFMM_LIB_DIR)/libpvfmm.a

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_BIN = \
       $(BINDIR)/dbl-test \
       $(BINDIR)/bie-solve

all : $(TARGET_BIN)

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $^ $(LDLIBS) -o $@ $(CXXFLAGS)
ifeq "$(OS)" "Darwin"
	/usr/bin/dsymutil $@ -o $@.dSYM
endif

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(SCTL_INCLUDE_DIR) -c $^ -o $@

.PHONY: all check clean

test: $(TARGET_BIN)
	for f in $^ ; do ./$$f ; done

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~ */*/*~
