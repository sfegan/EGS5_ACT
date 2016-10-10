// VSORayTracer.i - SWIG wrapper for the ChiLA (VERITAS) optics simulation
// Stephen Fegan - sfegan@llr.in2p3.fr - June 2013
// $Id: VSORayTracer.i 5681 2013-09-27 10:53:37Z sfegan $

%module VSORayTracer
%{
#include <iostream>
#include "Constants.hpp"
#include "Vec3D.hpp"
#include "Vec4D.hpp"
#include "Particle.hpp"
#include "VSOArrayParameters.hpp"
#include "VSOPixel.hpp"
#include "VSOMirror.hpp"
#include "VSOObscuration.hpp"
#include "VSORayTracer.hpp"
#include "VSOTelescope.hpp"
#include "VSOTelescopeArray.hpp"
#include "SimpleRNG.hpp"
#include "RandomNumbers.hpp"
#include "RandomNumbers_TNG.hpp"

using namespace std;
using namespace VERITAS;
using namespace Physics;
%}

%include "std_string.i"
%include "std_vector.i"
%include "Constants.hpp"
%include "Vec3D.hpp"
%include "Vec4D.hpp"
%include "Particle.hpp"

%{
  static Physics::Vec3D Vec3D_cross(const Physics::Vec3D& a, const Physics::Vec3D& b)
  {
    return a^b;
  }

  static double Vec3D_dot(const Physics::Vec3D& a, const Physics::Vec3D& b)
  {
    return a*b;
  }
  %}

Physics::Vec3D Vec3D_cross(const Physics::Vec3D& a, const Physics::Vec3D& b);
double Vec3D_dot(const Physics::Vec3D& a, const Physics::Vec3D& b);


%include "VSOArrayParameters.hpp"
%include "VSOMirror.hpp"
%include "VSOPixel.hpp"
%include "VSOObscuration.hpp"
%include "VSOTelescope.hpp"
%include "VSOTelescopeArray.hpp"
%include "VSORayTracer.hpp"
%include "SimpleRNG.hpp"
%include "RandomNumbers_TNG.hpp"
%template(RandomNumbers) RandomNumbers_TNG<RNGCore::NR3Ran>;
%template(SimpleRNGAdapter) SimpleRNGAdapter_TNG<RandomNumbers>;
%template(StdVectorDouble) std::vector<double>;
%include "RandomNumbers.hpp"
