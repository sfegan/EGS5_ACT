DEFINES   = -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS
OPT       = -g #-mssse3 -mfpmath=sse,387 -ftree-vectorize -march=core2
VERITASDIR = /Users/sfegan/Google\ Drive/Code/Projects/MyGlastCode/VERITAS
INCDIRS   =  -I$(VERITASDIR)
CFLAGS    = $(DEFINES) $(OPT) -Wall $(INCDIRS)
CXXFLAGS  = $(DEFINES) $(OPT) -Wall $(INCDIRS)
FFLAGS    = $(OPT) -fno-automatic -finit-local-zero -I$(VERITASDIR)
LDFLAGS   = $(OPT) -L/opt/local/lib -L$(VERITASDIR)
LIBS      = -lVERITAS -lhdf5 -lgfortran -lpthread

CC        = gcc-mp-4.6
CXX       = g++-mp-4.6
FC        = gfortran-mp-4.6
SWIG      = swig -cpperraswarn -python -c++ $(DEFINES) $(INCDIRS)
PYCONFIG  = python3.3-config

SWIGHEAD  = Atmosphere.hpp Cherenkov.hpp Interpolation1D.hpp \
		DetectorEfficiency.hpp BField.hpp EGS5System.hpp \
		EGS5LayeredDetector.hpp EGS5AtmosphericDetector.hpp \
		SimpleInstrumentedDetectors.hpp MultiRNG.hpp
FSRC      = egs5_system.f egs5_futils.f
CXXSRC    = Atmosphere.cpp EGS5System.cpp EGS5LayeredDetector.cpp \
		EGS5AtmosphericDetector.cpp BField.cpp DetectorEfficiency.cpp \
		SimpleInstrumentedDetectors.cpp MultiRNG.cpp
OBJ       = $(CXXSRC:.cpp=.o) $(FSRC:.f=.o)
PROGS     = image cherenkov_yield

export CFLAGS FFLAGS CXXFLAGS LDFLAGS CC CXX FC LIBS SWIG PYCONFIG

all: $(PROGS) lib

lib: _EGS5Simulations.so

test: test.o $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

cherenkov_yield: cherenkov_yield.o $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

image: image.o $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

mfp: mfp.o $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

egs5_system.f:
	cat egs5/egs/COPYRIGHT > egs5_system.f
	cat egs5/egs/*.f >> egs5_system.f
	cat egs5/auxcode/*.f >> egs5_system.f
	cat egs5/pegs/*.f >> egs5_system.f
	sed -i "" -e "s/include\//egs5\/include\//" egs5_system.f
	sed -i "" -e "s/auxcommons\//egs5\/auxcommons\//" egs5_system.f
	sed -i "" -e "s/pegscommons\//egs5\/pegscommons\//" egs5_system.f
	sed -i "" -e "s/data\//egs5\/data\//" egs5_system.f
	sed -i "" -e "s/CHARACTER\*24 FILE1,FILE2/CHARACTER\*29 FILE1,FILE2/" egs5_system.f
	sed -i "" -e "s/write(6,/write(66,/" egs5_system.f
	sed -i "" -e "s/subroutine randomset/subroutine origrandomset/" egs5_system.f
	awk 'BEGIN{c=0};{if(match($$0,/call randomset/)>0){split($$0,a,")");printf("%s,%i)%s\n",a[1],c,a[2]);c=c+1}else{print $$0;}};END{printf("      subroutine ncallrandomset(incall)\n      implicit none\n      integer incall\n      incall=%i\n      return\n      end\n",c)}' egs5_system.f > _egs5_system.f
	mv _egs5_system.f egs5_system.f

_EGS5Simulations.so: EGS5Simulations_wrap.o $(OBJ)
	$(CXX) `$(PYCONFIG) --ldflags` -shared $(LDFLAGS) -o $@ $^ $(LIBS) `$(PYCONFIG) --libs`

EGS5Simulations_wrap.o: EGS5Simulations_wrap.cxx
	$(CXX) $(CXXFLAGS) -I `$(PYCONFIG) --includes` -c $<

EGS5Simulations_wrap.cxx: EGS5Simulations.i $(SWIGHEAD)
	$(SWIG) $<

VERITAS/libVSOptics.a: RayTracing

SUBDIRS = RayTracing

$(SUBDIRS):
	$(MAKE) -C $@

$(addsuffix -clean,$(SUBDIRS)):
	$(MAKE) -C $(@:-clean=) clean

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

.PHONY: distclean clean $(SUBDIRS) $(addsuffix -clean,$(SUBDIRS))

clean: 
	$(RM) *.o $(PROGS) EGS5Simulations_wrap.cxx EGS5Simulations.pyEGS5Simulations.pyc _EGS5Simulations.so *~

distclean: clean $(addsuffix -clean,$(SUBDIRS))
	$(RM) egs5_system.f pgs5job.* egs5job.* scp.dat
