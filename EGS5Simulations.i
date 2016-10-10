// EGS5Simulations.i - SWIG wrapper for the EG5 simulation system
// Stephen Fegan - sfegan@llr.in2p3.fr - June 2013
// $Id: EGS5Simulations.i 5688 2013-09-30 13:58:43Z sfegan $

%module EGS5Simulations
%{
namespace VERITAS
{
};

#include <vector>
#include "VSSimple2DHist.hpp"
#include "VSAAlgebra.hpp"
#include "SimpleRNG.hpp"
#include "MultiRNG.hpp"
#include "Atmosphere.hpp"
#include "Cherenkov.hpp"
#include "DetectorEfficiency.hpp"
#include "BField.hpp"
#include "EGS5System.hpp"
#include "EGS5LayeredDetector.hpp"
#include "EGS5AtmosphericDetector.hpp"
#include "SimpleInstrumentedDetectors.hpp"
%}

namespace VERITAS
{
};

%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include "VSDataConverter.hpp"
%include "VSSimple2DHist.hpp"
%include "VSAAlgebra.hpp"
%include "Atmosphere.hpp"
%include "Cherenkov.hpp"
%include "Interpolation1D.hpp"
%include "SimpleRNG.hpp"
%include "MultiRNG.hpp"
%include "DetectorEfficiency.hpp"
%include "BField.hpp"
%include "EGS5System.hpp"
%template(VecMedia) std::vector<EGS5Media>;
%template(VecRegion) std::vector<EGS5Region>;
%include "EGS5LayeredDetector.hpp"
%template(VecLayer) std::vector<EGS5Layer>;
%include "EGS5AtmosphericDetector.hpp"
%template(VecScope) std::vector<EGS5ACTArrayScope>;
%template(VecImagingScope) std::vector<EGS5ACTArrayImagingScope>;
%template(Image) VERITAS::VSSimple2DHist<double, double>;
%include "SimpleInstrumentedDetectors.hpp"
%template(VecTrack) std::vector<Track>;
