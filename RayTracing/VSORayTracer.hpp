//-*-mode:c++; mode:font-lock;-*-

/*! \file VSORayTracer.hpp

  Ray tracing class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz             \n
           UCLA                        \n
	   nicewicz@physics.ucla.edu   \n

  \date    12/05/2004
  \version 0.2
  \note
*/

#ifndef VSORAYTRACER_HPP
#define VSORAYTRACER_HPP

#include <Vec3D.hpp>
#include <Vec4D.hpp>
#include <Particle.hpp>
#include <RandomNumbers.hpp>

#include "VSOTelescopeArray.hpp"

namespace VERITAS
{
  
  enum VSOTraceStatus { TS_NONE,                               // 0
			TS_DOES_INTERSECT_GROUND,              // 1
			TS_NO_SCOPE,                           // 2
			TS_OUTSIDE_REFLECTOR_IP,               // 3
			TS_MISSED_REFLECTOR_SPHERE,            // 4
			TS_NO_MIRROR,                          // 5
			TS_MIRROR_REMOVED,                     // 6
			TS_MISSED_MIRROR_SPHERE,               // 7
			TS_OBSCURED_BEFORE_MIRROR,             // 8
			TS_MISSED_MIRROR_EDGE,                 // 9
			TS_ABSORBED_AT_MIRROR,                 // 10
			TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE,   // 11
			TS_OBSCURED_BEFORE_FOCAL_PLANE,        // 12
			TS_NO_PIXEL,                           // 13
			TS_ABSORBED_AT_CONCENTRATOR,           // 14
			TS_PE_GENERATED };                     // 15

  class VSOTraceInfo
  {
  public:					
    const VSOTelescopeArray* array;
    VSOTraceStatus      status;
    double              ground_x;
    double              ground_y;
    double              ground_dx;
    double              ground_dy;
    int                 scope_hexid;
    const VSOTelescope* scope;
    double              reflec_x;
    double              reflec_z;
    double              hex_reflec_x;
    double              hex_reflec_z;
    double              hex_reflec_dx;
    double              hex_reflec_dz;
    int                 mirror_hexid_nominal;
    int                 mirror_hexid;
    const VSOMirror*    mirror;
    double              mirror_x;
    double              mirror_y;
    double              mirror_z;
    Physics::Vec3D      mirror_normal;
    double              mirror_normal_dispersion;
    Physics::Vec3D      mirror_scattered;
    double              mirror_reflection_angle;
    double              fplane_x;
    double              fplane_z;
    double              fplane_dx;
    double              fplane_dz;
    double              fplane_t;
    int                 pixel_hexid;
    const VSOPixel*     pixel;
    double              pixel_dist;
    bool                concentrator_hit;
    unsigned            obscuration_id;
    const VSOObscuration* obscuration;
    
    void reset();
    std::ostream& write(std::ostream& stream = std::cout, bool convert_to_physical_units=false) const;
    
    bool rayWasReflected() const 
    { return (int)status >= (int)TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE; }
    bool rayHitFocalPlane() const 
    { return (int)status >= (int)TS_NO_PIXEL; }
    
  };

  class VSOPSFInfo
  {
  public:					
    unsigned            nhit;
    std::vector<double> r_tan;
    std::vector<double> r_sag;
    std::vector<double> t;
    double              mean_tan;
    double              mean_sag;
    double              mean_t;
    double              rms_tan;
    double              rms_sag;
    double              cov_tan_sag;
    double              rms_t;
    double              median_tan;
    double              median_sag;
    double              median_t;
    double              r80;
    double              t80;

    void reset() { *this = VSOPSFInfo(); }
  };

  class VSORayTracer
  {
  public:
    VSORayTracer(const VSOTelescopeArray& array, RandomNumbers& rng): 
      fArray(array), fRNG(rng) { }
    virtual ~VSORayTracer();
    
    typedef VSOTraceStatus Status;
    typedef VSOTraceInfo TraceInfo;
    
    const VSOPixel* trace(Physics::Particle& ray, TraceInfo& info);
    const VSOPixel* trace(Physics::Particle& ray, TraceInfo& info,
			  const VSOTelescope* scope_hint);
    
    bool beam(Physics::Particle& photon,
	      const Physics::Vec3D& origin, const Physics::Vec3D& direction, 
	      double beam_start, double beam_stop, 
	      double beam_radius_in, double beam_radius_out,
	      double beam_angle_lo, double beam_angle_hi,
	      double lambda_nm = 400);
    
    bool laserBeam(Physics::Particle& photon, 
		   const Physics::Vec3D& center,
		   const Physics::Vec3D& direction, 
		   double d0, double sampling_radius,
		   double lambda_nm = 400);
    
    bool fanBeam(Physics::Particle& photon,
		 const Physics::Vec3D& origin, 
		 const Physics::Vec3D& direction, 
		 double half_angle_spread, double lambda_nm = 400);
    
    bool muonBeam(Physics::Particle& photon,
		  const Physics::Vec3D& origin, 
		  const Physics::Vec3D& direction, 
		  double muon_travel_distance, double opening_angle, 
		  double lambda_nm = 400);

    bool testBeam(Physics::Particle& photon,
		  const VSOTelescope* scope,
		  double theta, double phi = 0,
		  double U = std::numeric_limits<double>::infinity(),
		  double lambda_nm = 400);

    void calcPSF(class VSOPSFInfo& psf, const VSOTelescope* scope,
		 double theta, double phi = 0,
		 double U = std::numeric_limits<double>::infinity(),
		 unsigned nsim = 1000000, bool save_image = false);

  private:
    bool findMirror(Physics::Particle& ray, TraceInfo& info);

    const VSOTelescopeArray& fArray;
    RandomNumbers&           fRNG;

    const VSOPixel* scope_trace(Physics::Particle& ray, TraceInfo& info);
  };
  
#ifndef SWIG
  std::ostream& operator <<(std::ostream& stream, 
			    const VSORayTracer::TraceInfo& o);
#endif
  
}

#endif // VSORAYTRACER_HPP
