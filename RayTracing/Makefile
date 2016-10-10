LIBOBJECTS = VSOArrayParameters.o VSOObscuration.o \
	     VSOMirror.o VSOPixel.o VSOTelescope.o VSOTelescopeArray.o \
             VSORayTracer.o Vec3D.o Vec4D.o Particle.o Constants.o xytohex.o

BINOBJECTS = test_raytrace.o

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVSOptics.a _VSORayTracer.so $(BINOBJECTS:.o=)

CXXFLAGS += -I.
LDFLAGS += -L.
LIBS := -lVSOptics $(LIBS)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSOptics.a: $(LIBOBJECTS)
	$(AR) r $@ $^

test_raytrace: test_raytrace.o libVSOptics.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

VSOArrayParameters.o: VSOArrayParameters.cpp VSOArrayParameters.hpp

VSOTelescopeArray.o: VSOArrayParameters.hpp

VSOArrayParameters.cpp: ARRAY_INI_TEMPLATE code_gen_array_ini_cpp.pl
	perl code_gen_array_ini_cpp.pl $< > $@

VSOArrayParameters.hpp: ARRAY_INI_TEMPLATE code_gen_array_ini_hpp.pl
	perl code_gen_array_ini_hpp.pl $< > $@

_VSORayTracer.so: VSORayTracer_wrap.o $(LIBOBJECTS)
	$(CXX) `$(PYCONFIG) --ldflags` $(LDFLAGS) -shared -o $@ $^ $(LIBS) `$(PYCONFIG) --libs`

VSORayTracer_wrap.o: VSORayTracer_wrap.cxx
	$(CXX) $(CXXFLAGS) -I `$(PYCONFIG) --includes` -c $<

VSORayTracer_wrap.cxx: VSORayTracer.i $(LIBOBJECTS:.o=.hpp)
	$(SWIG) $<

.PHONY: clean

clean:
	$(RM) \
	VSOArrayParameters.cpp VSOArrayParameters.hpp \
	_VSORayTracer.so VSORayTracer_wrap.o VSORayTracer_wrap.cxx \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
