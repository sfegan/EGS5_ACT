//-*-mode:c++; mode:font-lock;-*-

/*! \file RayTracer.cpp

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

#include <cmath>
#include <iostream>
#include <algorithm>

#include <Vec3D.hpp>
#include <Vec4D.hpp>

#include <VSDataConverter.hpp>

#include "VSORayTracer.hpp"

using namespace Physics;
using namespace VERITAS;

void VSOTraceInfo::reset()
{
  *this = VSOTraceInfo();
}

std::ostream& VSOTraceInfo::write(std::ostream& stream, bool cpu) const
{
  stream 
    << VSDataConverter::toString(status) << ' '                                                    // $1  --  AWK column
    << VSDataConverter::toString(ground_x * ( cpu ? array->spacing() : 1.0 )) << ' '               // $2
    << VSDataConverter::toString(ground_y * ( cpu ? array->spacing() : 1.0 )) << ' '               // $3
    << VSDataConverter::toString(ground_dx * ( cpu ? array->spacing() : 1.0 )) << ' '              // $4
    << VSDataConverter::toString(ground_dy * ( cpu ? array->spacing() : 1.0 )) << ' '              // $5
    << VSDataConverter::toString(scope_hexid) << ' '                                               // $6
    << VSDataConverter::toString(( scope ? scope->hexID() : 0 )) << ' '                            // $7
    << VSDataConverter::toString(reflec_x * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $8
    << VSDataConverter::toString(reflec_z * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $9
    << VSDataConverter::toString(hex_reflec_x * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $10
    << VSDataConverter::toString(hex_reflec_z * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $11
    << VSDataConverter::toString(hex_reflec_dx * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '  // $12
    << VSDataConverter::toString(hex_reflec_dz * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '  // $13
    << VSDataConverter::toString(mirror_hexid) << ' '                                              // $14
    << VSDataConverter::toString(( mirror ? mirror->hexID() : 0 )) << ' '                          // $15
    << VSDataConverter::toString(mirror_x) << ' '                                                  // $16
    << VSDataConverter::toString(mirror_y) << ' '                                                  // $17
    << VSDataConverter::toString(mirror_z) << ' '                                                  // $18
    << VSDataConverter::toString(mirror_normal) << ' '                                             // $19 $20 $21 $22 $23
    << VSDataConverter::toString(mirror_normal_dispersion) << ' '                                  // $24
    << VSDataConverter::toString(mirror_scattered) << ' '                                          // $25 $26 $27 $28 $29
    << VSDataConverter::toString(mirror_reflection_angle) << ' '                                   // $30
    << VSDataConverter::toString(fplane_x * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '   // $31
    << VSDataConverter::toString(fplane_z * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '   // $32
    << VSDataConverter::toString(fplane_dx * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '  // $33
    << VSDataConverter::toString(fplane_dz * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '  // $34
    << VSDataConverter::toString(fplane_t) << ' '                                                  // $35
    << VSDataConverter::toString(pixel_hexid) << ' '                                               // $36
    << VSDataConverter::toString(( pixel ? pixel->hexID() : 0 )) << ' '                            // $37
    << VSDataConverter::toString(pixel_dist * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' ' // $38
    << VSDataConverter::toString(concentrator_hit);                                                // $39
  return stream;
}


std::ostream& operator << (std::ostream& stream, const VSOTraceInfo& o)
{
  return o.write(stream,false);
}


VSORayTracer::~VSORayTracer()
{
  // nothing to see here
}

const VSOPixel* VSORayTracer::trace(Physics::Particle& ray, TraceInfo& info,
				    const VSOTelescope* scope_hint)
{
  // Initialize array
  info.reset();
  info.array = &fArray;

  info.scope = scope_hint;

  return(scope_trace(ray,info));
}

const VSOPixel* VSORayTracer::trace(Physics::Particle& ray, TraceInfo& info)
{
  // Initialize array
  info.reset();
  info.array = &fArray;

  // Require photon to be down going
  if(ray.Momenta().r.z>0)
    {
      info.status = TS_DOES_INTERSECT_GROUND;
      return 0;
    }

  // Propagate to ground
  Particle ray_copy(ray);
  bool good;
  good = ray_copy.PropagateFreeToPlane(Vec3D(0,0,1), -fArray.altitude(), true);
  if(!good)
    {
      info.status = TS_DOES_INTERSECT_GROUND;
      return 0;
    }

  // Array is assumed to be hexagonal use VVV look-up function to find site
  info.ground_x = ray_copy.Position().r.x / fArray.spacing();
  info.ground_y = ray_copy.Position().r.y / fArray.spacing();
  info.ground_dx = info.ground_x;
  info.ground_dy = info.ground_y;
  if(fArray.arrayParity())info.ground_dx = -info.ground_dx;
  xy_to_nh(&info.ground_dx,&info.ground_dy,&info.scope_hexid);
  if(fArray.arrayParity())info.ground_dx = -info.ground_dx;

  // Find telescope (if there is a real telescope at that site)
  info.scope = fArray.telescopeByHexID(info.scope_hexid);
  if(info.scope==0)
    {
      info.status = TS_NO_SCOPE;
      return 0;
    }

#if 0 
  // Propagate to reflector impact sphere -- test whether ray can hit scope
  // This probably does not speed things up much as the propagation to
  // the reflector (two steps down) does a similar thing
  good = ray_copy.PropagateFreeToSphere(info.scope->position(),
					info.scope->reflectorIP(),
					IP_EARLIEST,true);
  if(!good)
    {
      info.status = TS_OUTSIDE_REFLECTOR_IP;
      return 0;
    }
#endif

  return(scope_trace(ray,info));
}

const VSOPixel* 
VSORayTracer::scope_trace(Physics::Particle& ray, TraceInfo& info)
{
  bool good;

  // Transform to reflector co-ordinate system
  info.scope->globalToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "A: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif
  // **************************************************************************
  // ****************** RAY IS NOW IN RELECTOR COORDINATES *******************
  // **************************************************************************

  // Test for obscuration
  unsigned nobs = info.scope->numObscurations();
  unsigned obs_ihit  = nobs;
  double   obs_time  = ray.Position().r0;
  for(unsigned iobs=0;iobs<nobs;iobs++)
    {
      Particle p_out;
      if(info.scope->obscuration(iobs)->doesObscure(ray, p_out))
	{
	  if((obs_ihit==nobs)||(p_out.Position().r0<obs_time))
	    obs_ihit = iobs, obs_time = p_out.Position().r0;
	}
    }

  // Propagate to intersection with the reflector sphere
  good = ray.PropagateFreeToSphere(Vec3D(0,info.scope->curvatureRadius(),0),
				   info.scope->curvatureRadius(), 
				   Particle::IP_LATEST, false /* true */);
  if(!good)
    {
      info.status = TS_MISSED_REFLECTOR_SPHERE;
      return 0;
    }

  // Assume mirrors on hexagonal grid - use VVV routines to find which mirror
  double tx = ray.Position().r.x / info.scope->facetSpacing();
  double tz = ray.Position().r.z / info.scope->facetSpacing();
  info.reflec_x     = tx;
  info.reflec_z     = tz;
  double costheta = cos(info.scope->reflectorRotation());
  double sintheta = sin(info.scope->reflectorRotation());
  info.hex_reflec_x = tx*costheta - tz*sintheta;  // Align with hex grid in dir of
  info.hex_reflec_z = tz*costheta + tx*sintheta;  // Vec3D(0,-reflectorRotation,0)
  info.hex_reflec_dx = info.hex_reflec_x;
  info.hex_reflec_dz = info.hex_reflec_z;
  if(info.scope->mirrorParity())info.hex_reflec_dx = -info.hex_reflec_dx; // Reverse parity if required
  xy_to_nh(&info.hex_reflec_dx,&info.hex_reflec_dz,&info.mirror_hexid);
  info.mirror_hexid_nominal = info.mirror_hexid;
  if(info.scope->mirrorParity())info.hex_reflec_dx = -info.hex_reflec_dx; // Reverse parity if required
  
  // Find mirror - searching neighbors if desired  
  good = findMirror(ray, info);
  if(!good)
    {
      // info.status should be set appropriately by findMirror
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // **************************************************************************
  // ******************** RAY IS NOW IN MIRROR COORDINATES ********************
  // **************************************************************************

  // Test for interaction with obscuration before mirror was hit
  if((obs_ihit!=nobs) && (obs_time<ray.Position().r0))
    {
      info.status = TS_OBSCURED_BEFORE_MIRROR;
      info.mirror->mirrorToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      info.obscuration = info.scope->obscuration(obs_ihit);
      info.obscuration_id = obs_ihit;
      return 0;
    }

  // Check to see if photon is absorbed at the mirror. Would be faster to 
  // do this check before the edge check, but this way gets info.status 
  // correct .. is this important ?
  if(fRNG.Uniform() > info.mirror->degradingFactor())
    {
      info.status = TS_ABSORBED_AT_MIRROR;
      info.mirror->mirrorToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Find normal at the point of intersection
  double mirror_radius = info.mirror->focalLength()*2.0;
  info.mirror_normal = Vec3D(0,mirror_radius,0) - ray.Position().r;
  info.mirror_normal /= info.mirror_normal.Norm();

  // Scatter the normal to account for the spot size ot the focal length of the
  // radius. The spot size is given as the DIAMETER at the focal distance.
  // Must divide by 2.0 (for reflection)
  info.mirror_normal_dispersion = 
    info.mirror->spotSize()/2.0/info.mirror->focalLength();

  info.mirror_scattered = info.mirror_normal;
  info.mirror_scattered.ScatterDirection(info.mirror_normal_dispersion,fRNG);
 
  // Reflect ray
  ray.Reflect(info.mirror_scattered);
  
  // Back to reflector coordinates
  info.mirror->mirrorToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "C: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif

  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  // Test for obscuration
  obs_ihit  = nobs;
  obs_time  = ray.Position().r0;
  for(unsigned iobs=0;iobs<nobs;iobs++)
    {
      Particle p_out;
      if(info.scope->obscuration(iobs)->doesObscure(ray, p_out))
	{
	  if((obs_ihit==nobs)||(p_out.Position().r0<obs_time))
	    obs_ihit = iobs, obs_time = p_out.Position().r0;
	}
    }

  // Translate to focal plane coordinates
  info.scope->reflectorToFocalPlane(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "D: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif

  // **************************************************************************
  // ***************** RAY IS NOW IN FOCAL PLANE COORDINATES ******************
  // **************************************************************************

  // Propagate back to camera plane
  good = ray.PropagateFreeToPlane(Vec3D(0,1,0),0,false);
  if(!good)
    {
      info.status = TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE;
      info.scope->focalPlaneToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Test for interaction with obscuration before focal plane was hit
  if((obs_ihit!=nobs) && (obs_time<ray.Position().r0))
    {
      info.status = TS_OBSCURED_BEFORE_FOCAL_PLANE;
      info.scope->focalPlaneToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      info.obscuration = info.scope->obscuration(obs_ihit);
      info.obscuration_id = obs_ihit;
      return 0;
    }

  info.fplane_x = ray.Position().r.x / info.scope->pixelSpacing();
  info.fplane_z = ray.Position().r.z / info.scope->pixelSpacing();
  info.fplane_dx = info.fplane_x;
  info.fplane_dz = info.fplane_z;
  info.fplane_t = ray.Position().r0 / Constants::cgs_c;
  if(info.scope->pixelParity())info.fplane_dx = -info.fplane_dx;
  xy_to_nh(&info.fplane_dx,&info.fplane_dz,&info.pixel_hexid);
  if(info.scope->pixelParity())info.fplane_dx = -info.fplane_dx;

  // Find pixel (if there is a real pixel at that site)
  info.pixel = info.scope->pixelByHexID(info.pixel_hexid);
  if(info.pixel==0)
    {
      info.status = TS_NO_PIXEL;
      info.scope->focalPlaneToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  info.pixel_dist = 
    sqrt(info.fplane_dx*info.fplane_dx + info.fplane_dz*info.fplane_dz) * 
    info.scope->pixelSpacing();

  info.concentrator_hit = 
    (info.pixel_dist > info.scope->cathodeDiameter()/2.0);
  if(info.concentrator_hit)
    {
      if(fRNG.Uniform() > info.scope->concentratorSurvivalProb())
	{
	  info.status = TS_ABSORBED_AT_CONCENTRATOR;
	  info.scope->focalPlaneToReflector(ray);
	  info.scope->reflectorToGlobal(ray);
	  return 0;
	}
    }

  // Translate to reflector coordinates
  info.scope->focalPlaneToReflector(ray);

  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  info.pixel_dist = Vec3D(ray.Position().r - info.pixel->pos() -
			  info.scope->focalPlanePosition()).Norm();

  // Transform back to global
  info.scope->reflectorToGlobal(ray);

  // **************************************************************************
  // ******************** RAY IS NOW IN GLOBAL COORDINATES ********************
  // **************************************************************************
    
  info.status = TS_PE_GENERATED;
  return info.pixel;
}

// Find the mirror (allowing that it might not actually be one that
// neighbors the optimal position). Starts with ray at the reflector
// sphere, in the reflector coordinate system and info.mirror_hexid
// set to "nominal" mirror given position on reflector. Ends with
// either the ray in the mirror system (with info.mirror,
// mirror_hexid, mirror_x, mirror_y, mirror_z set correctly), in which
// case "true" is returned, or the mirror still in reflector
// coordinates, in the case where no mirror is encountered, where
// false is returned

bool VSORayTracer::findMirror(Physics::Particle& ray, TraceInfo& info)
{
  // **************************************************************************
  // ****************** RAY STARTS IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  const static unsigned search_radius = 0;
  const static unsigned nsearch = 3*search_radius*(search_radius+1)+1;

  Physics::Particle ray_in(ray);

#if 0
  int               info_hexid  = info.mirror_hexid;
  Status            info_status = info.status;
  const VSOMirror*  info_mirror;
  bool              info_mirror_x;
  bool              info_mirror_y;
  bool              info_mirror_z;
  Particle          info_ray;
#endif
  bool              impinging_ray_found = false;

  int hex_i0 = info.mirror_hexid;
  double hex_x0 = 0;
  double hex_y0 = 0;
  if(nsearch > 1)nh_to_xy(&hex_i0, &hex_x0, &hex_y0);
  
  for(unsigned isearch = 0; isearch<nsearch; isearch++)
    {
      int test_hexid;
      if(isearch == 0)
	{
	  test_hexid = info.mirror_hexid;
	}
      else
	{
	  int hex_i1 = isearch;
	  double hex_x1 = 0;
	  double hex_y1 = 0;
	  nh_to_xy(&hex_i1, &hex_x1, &hex_y1);
	  hex_x1 += hex_x0;
	  hex_y1 += hex_y0;
	  xy_to_nh(&hex_x1,&hex_x0,&test_hexid);
	}

      // Find mirror (if there is a real mirror at that site)
      const VSOMirror* test_mirror = info.scope->mirrorByHexID(test_hexid);
      if(isearch==0)info.mirror = test_mirror;

      if(test_mirror==0)
	{
	  if(isearch==0)info.status = TS_MIRROR_REMOVED;
	  continue;
	}

      if(test_mirror->removed())
	{
	  if(isearch==0)info.status = TS_MIRROR_REMOVED;
	  continue;
	}

      // Propagate to intersection with the mirror sphere
      double mirror_radius = test_mirror->focalLength()*2.0;
      Vec3D mirror_center = 
	test_mirror->pos() + test_mirror->align()*mirror_radius;

      Particle test_ray(ray_in);
      bool good;
      good = test_ray.PropagateFreeToSphere(mirror_center, mirror_radius,
					    Particle::IP_LATEST, true);
      if(isearch==0)ray = test_ray;
      if(!good)
	{
	  if(isearch==0)info.status = TS_MISSED_MIRROR_SPHERE;
	  continue;
	}

      // Convert to mirror coordinates
      test_mirror->reflectorToMirror(test_ray);
      if(isearch==0)ray = test_ray;

      // **********************************************************************
      // ****************** RAY IS NOW IN MIRROR COORDINATES ******************
      // **********************************************************************

      double test_mirror_x = test_ray.Position().r.x;
      double test_mirror_y = test_ray.Position().r.y;
      double test_mirror_z = test_ray.Position().r.z;

      if(isearch==0)
	{
	  info.mirror_x  = test_mirror_x;
	  info.mirror_y  = test_mirror_y;
	  info.mirror_z  = test_mirror_z;
	}
      
      // Check if ray impacted beyond the edge of this mirror
      static const double cos60 = 1.0/2.0;
      static const double sin60 = sqrt(3.0)/2.0;
      double edge = info.scope->facetSize()/2.0;
      double x_0 = test_mirror_x;
      double x_pos60 = cos60*test_mirror_x - sin60*test_mirror_z;
      double x_neg60 = cos60*test_mirror_x + sin60*test_mirror_z;

      if((x_0>edge)||(x_0<-edge)||(x_pos60>edge)||(x_pos60<-edge)||
	 (x_neg60>edge)||(x_neg60<-edge))
	{
	  if(isearch==0)info.status = TS_MISSED_MIRROR_EDGE;
	  continue;
	}
  
      // We have a good ray that impinges on a mirror - test that it hits
      // earlier than previously found one (if any)
      if((isearch>0)
	 &&((!impinging_ray_found)
	    ||(test_ray.Position().r0<ray.Position().r0)))
	{
	  info.mirror_hexid = test_hexid;
	  info.status       = info.status;
	  info.mirror       = test_mirror;
	  info.mirror_x     = test_mirror_x;
	  info.mirror_y     = test_mirror_y;
	  info.mirror_z     = test_mirror_z;
	  ray               = test_ray;
	}
      impinging_ray_found = true;
    }

  return impinging_ray_found;
}

bool VSORayTracer::beam(Physics::Particle& photon,
		     const Physics::Vec3D& origin,
		     const Physics::Vec3D& direction, 
		     double beam_start, double beam_stop, 
		     double beam_radius_in, double beam_radius_out,
		     double beam_angle_lo, double beam_angle_hi,
		     double lambda_nm)
{
  double d_norm = direction.Norm();
  if(d_norm == 0)return false;
  Physics::Vec3D d_hat = direction/d_norm;

  Vec3D tangent_a(1,0,0);
  Vec3D tangent_b(0,0,1);

  Vec3D rot = Vec3D(0,-1,0)^d_hat;
  if(rot.Norm())
    {
      rot *= atan2(rot.Norm(), Vec3D(0,-1,0)*d_hat)/rot.Norm();
      tangent_a.Rotate(rot);
      tangent_b.Rotate(rot);
    }

  // CHOOSE PHOTON EMISSION POINT
  Physics::Vec3D emission_point(origin);
  emission_point += d_hat*beam_start;

  // SAMPLE PHOTON EMISSION POINT FROM BEAM LENGTH
  if(beam_stop != beam_start)
    {
      emission_point += d_hat*((beam_stop-beam_start)*fRNG.Uniform());
    }

  // SAMPLE PHOTON EMISSION POINT FROM BEAM AREA
  if((beam_radius_in != 0)||(beam_radius_out != 0))
    {
      double theta = 2.0*M_PI*fRNG.Uniform();
      double rho2 = beam_radius_in*beam_radius_in;
      if(beam_radius_in != beam_radius_out)
	rho2 += (beam_radius_out*beam_radius_out-rho2) * fRNG.Uniform();
      double rho = sqrt(rho2);
      emission_point += tangent_a*rho*cos(theta) + tangent_b*rho*sin(theta);
    }

  // SAMPLE PHOTON EMISSION DIRECTION
  double costheta = cos(beam_angle_lo);
  if(beam_angle_lo != beam_angle_hi)
    {
      costheta += fRNG.Uniform() * (cos(beam_angle_hi)-costheta);
    }
  double theta = acos(costheta);

  if(theta != 0)
    {
      double phi = fRNG.Uniform() * 2.0*Constants::num_pi;
#if 1
      d_hat = d_hat*costheta + (tangent_a*cos(phi)+tangent_b*sin(phi))*sin(theta);
#else  
      Vec3D axis = tangent_a*cos(phi) + tangent_b*sin(phi);
      d_hat.Rotate(axis*theta);
#endif
    }

  // MAKE PHOTON
  double E = Constants::cgs_h * Constants::cgs_c / ( lambda_nm * 1e-7 );
  photon = Particle(Physics::Vec4D(0,emission_point),d_hat,E,0,0);
  return true;
}

bool VSORayTracer::laserBeam(Physics::Particle& photon,
			  const Physics::Vec3D& origin,
			  const Physics::Vec3D& direction,
			  double d0, double sampling_radius,
			  double lambda_nm)
{
  return beam(photon, origin, direction, d0, d0, 0, sampling_radius, 0, 0,
	      lambda_nm);
}

bool VSORayTracer::fanBeam(Physics::Particle& photon,
			const Physics::Vec3D& origin, 
			const Physics::Vec3D& direction, 
			double half_angle_spread, 
			double lambda_nm)
{
  return beam(photon, origin, direction, 0, 0, 0, 0, 0, half_angle_spread,
	      lambda_nm);
}

bool VSORayTracer::muonBeam(Physics::Particle& photon,
			 const Physics::Vec3D& origin, 
			 const Physics::Vec3D& direction, 
			 double muon_travel_distance, double opening_angle, 
			 double lambda_nm)
{
  return beam(photon, origin, direction, 0, muon_travel_distance, 0, 0, 
	      opening_angle, opening_angle, lambda_nm);
}

bool VSORayTracer::testBeam(Physics::Particle& photon,
			    const VSOTelescope* scope,
			    double theta, double phi,
			    double U, double lambda_nm)
{
  if(std::isinf(U))
    {
      // Parallel beam - must set sampling radius
      Vec3D beam_dir(-sin(theta)*sin(phi) ,-cos(theta), -sin(theta)*cos(phi));
      Vec3D beam_cen(0, 0, 0);
      scope->reflectorToGlobal_mom(beam_dir);
      scope->reflectorToGlobal_pos(beam_cen);
      return laserBeam(photon, beam_cen, beam_dir, -
		       2.0*scope->focalPlanePosition().y, 
		       0.5*scope->reflectorIP(), lambda_nm);
    }
  else
    {
      double Utantheta = U*tan(theta);
      Vec3D beam_cen(Utantheta*sin(phi), U, Utantheta*cos(phi));
      double dist = beam_cen.Norm();
      Vec3D beam_dir(beam_cen*(-1.0/dist));
      scope->reflectorToGlobal_mom(beam_dir);
      scope->reflectorToGlobal_pos(beam_cen);
      return fanBeam(photon, beam_cen, beam_dir, 
		     asin(0.5*scope->reflectorIP()/dist), lambda_nm);
    }
}

void VSORayTracer::calcPSF(class VSOPSFInfo& psf, const VSOTelescope* scope,
			   double theta, double phi, double U, unsigned nsim,
			   bool save_image)
{
  double PS = scope->pixelSpacing()/scope->focalPlanePosition().y/M_PI*180.0;
  psf.reset();
  Particle ph;
  VSOTraceInfo info;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> t;
  std::vector<double> r2;  
  x.reserve(nsim);
  y.reserve(nsim);
  t.reserve(nsim);
  double cphi = std::cos(phi);
  double sphi = std::sin(phi);
  unsigned nhit = 0;
  double sumx   = 0;
  double sumy   = 0;
  double sumt   = 0;
  double sumxx  = 0;
  double sumyy  = 0;
  double sumxy  = 0;
  double sumtt  = 0;
  for(unsigned isim=0;isim<nsim;isim++)
    {
      testBeam(ph, scope, theta, phi, U);
      trace(ph, info, scope);
      if(!info.rayHitFocalPlane())continue;
      double _x = -PS*(cphi*info.fplane_z+sphi*info.fplane_x);
      double _y = -PS*(cphi*info.fplane_x-sphi*info.fplane_z);
      double _t = info.fplane_t;
      double x2 = _x*_x;
      double y2 = _y*_y;
      nhit++;
      sumx  += _x;
      sumy  += _y;
      sumt  += _t;
      sumxx += x2;
      sumxy += _x*_y;
      sumyy += y2;
      sumtt += _t*_t;
      x.push_back(_x);
      y.push_back(_y);
      t.push_back(_t);
    }

  psf.nhit        = nhit;

  if(save_image)
    {
      psf.r_tan   = x;
      psf.r_sag   = y;
      psf.t       = t;
    }

  psf.mean_tan    = sumx/double(nhit);
  psf.mean_sag    = sumy/double(nhit);
  psf.mean_t      = sumt/double(nhit);
  psf.rms_tan     = std::sqrt(sumxx/double(nhit) - psf.mean_tan*psf.mean_tan);
  psf.rms_sag     = std::sqrt(sumyy/double(nhit) - psf.mean_sag*psf.mean_sag);
  psf.cov_tan_sag = sumxy/double(nhit) - psf.mean_tan*psf.mean_sag;
  psf.rms_t       = std::sqrt(sumtt/double(nhit) - psf.mean_t*psf.mean_t);

  r2.resize(x.size());
  for(unsigned i=0;i<x.size(); i++)
    {
      double _x = x[i] - psf.mean_tan;
      double _y = y[i] - psf.mean_sag;
      r2[i] = _x*_x + _y*_y;
    }

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  std::sort(t.begin(), t.end());
  
  unsigned imed = x.size()/2;
  psf.median_tan  = x[imed];
  psf.median_sag  = y[imed];
  psf.median_t    = t[imed];

  for(unsigned i=0;i<x.size(); i++)
    {
      t[i] = std::fabs(t[i] - psf.mean_t);
    }

  std::sort(t.begin(), t.end());
  std::sort(r2.begin(), r2.end());

  unsigned i80 = (x.size()*8)/10;
  psf.r80         = std::sqrt(r2[i80]);
  psf.t80         = t[i80];
}
