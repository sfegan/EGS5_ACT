// EGS5AtmosphericDetector.hpp - Detector class for EGS5 system that
// - implements a layered atmospheric detector
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2012
// $Id: EGS5AtmosphericDetector.cpp 4954 2013-01-18 08:44:50Z sfegan $

#include <fstream>
#include <iostream>
#include <iomanip>

#include "EGS5AtmosphericDetector.hpp"

// ----------------------------------------------------------------------------
//
// EGS5AtmosphericDetector
//
// ----------------------------------------------------------------------------

EGS5AtmosphericDetector::
EGS5AtmosphericDetector(Atmosphere& atm, unsigned nlayer, double emax,
			BField* bfield, double zbottom, double ztop,
			unsigned nmedia)
: EGS5LayeredDetector(), m_atm(), m_bfield(), m_air_media()
{
  initialize(&atm, nlayer, emax, bfield, zbottom, ztop, nmedia);
}

EGS5AtmosphericDetector::~EGS5AtmosphericDetector()
{
  // nothing to see here
}

void EGS5AtmosphericDetector::
initialize(Atmosphere* atm, unsigned nlayer, double emax,
	   BField* bfield, double zbottom, double ztop, unsigned nmedia)
{
  m_atm = atm;
  m_bfield = bfield;
  std::vector<Media> media;
  std::vector<Layer> layers;
  makeMediaAndLayers(media, layers, *atm, nlayer, bfield, zbottom, ztop,
		     nmedia);
  EGS5LayeredDetector::initialize(media, layers, zbottom, emax);

  m_air_media.resize(media.size());
  for(unsigned imedia=0;imedia<media.size();imedia++)
    {
      m_air_media[imedia].first = media[imedia];
      m_air_media[imedia].second = 0;
    };
  for(unsigned ilayer=0;ilayer<layers.size();ilayer++)
    {
      unsigned imedia = layers[ilayer].imedia;
      m_air_media[imedia].second = 
	std::max(m_air_media[imedia].second, layers[ilayer].rho);
    }
}

void EGS5AtmosphericDetector::
makeMediaAndLayers(std::vector<Media>& media, std::vector<Layer>& layers,
		   Atmosphere& atm, unsigned nlayer, 
		   BField* bfield, double zbottom, double ztop,
		   unsigned nmedia)
{
  std::vector<AtmSlice> atm_slices;
  if(ztop == HUGE_VAL)ztop = atm.topOfAtmosphere();
  atm.makeAtmSlices(atm_slices, nlayer, ztop, zbottom);

  media.clear();

  layers.resize(atm_slices.size());
  bool has_bfield = (bfield!=0);
  double zb = zbottom;
  for(unsigned ilayer=0;ilayer<atm_slices.size();ilayer++)
    {
      double bx = 0;
      double by = 0;
      double bz = 0;
      double zt = atm_slices[ilayer].zt;

      unsigned imedia = (ilayer*nmedia)/atm_slices.size();

      if(has_bfield)bfield->getFieldCGS(0.5*(zb+zt), bx, by, bz);      
      layers[ilayer].imedia = imedia;
      layers[ilayer].zt     = zt;
      layers[ilayer].rho    = atm_slices[ilayer].rho;
      layers[ilayer].bfield = has_bfield;
      layers[ilayer].bx     = bx;
      layers[ilayer].by     = by;
      layers[ilayer].bz     = bz;
      zb = zt;

      if(imedia >= media.size())
	{
	  std::ostringstream name;
	  name << "AIR MEDIUM " << imedia+1;
	  media.push_back(Media(name.str(),1000));
	}
    }  
}

void EGS5AtmosphericDetector::
writePEGS5InputFile(const std::string& filename)
{
  std::ofstream stream(filename.c_str());
  for(unsigned imedia=0; imedia<m_air_media.size(); imedia++)
    stream
      << "MIXT\n"
      << " &INP NE=4,RHO=" << m_air_media[imedia].second
      << ",RHOZ=0.755215,0.231772,0.012882,0.000131\n"
      << "      IRAYL=0,IBOUND=0,INCOH=0,ICPROF=0,IMPACT=0 &END\n"
      << std::left << std::setw(30) << m_air_media[imedia].first.name
      << "AIR-GAS          \n"
      << "N  O  AR C \n"
      << "ENER\n"
      << " &INP AE=20.0,AP=20.0,UE=1.0E8,UP=1.0E8 &END\n"
      << "PWLF\n"
      << " &INP  &END\n"
      << "DECK\n"
      << " &INP  &END\n";
}

// ----------------------------------------------------------------------------
//
// EGS5AtmosphericCherenkovDetector
//
// ----------------------------------------------------------------------------

EGS5AtmosphericCherenkovDetector::
EGS5AtmosphericCherenkovDetector(Atmosphere& atm, unsigned nlayer, 
				 double emax, BField* bfield, 
				 double zbottom, double ztop, unsigned nmedia):
EGS5AtmosphericDetector(atm, nlayer, emax, bfield, zbottom, ztop, nmedia)
{
  // nothing to see here
}

EGS5AtmosphericCherenkovDetector::~EGS5AtmosphericCherenkovDetector()
{
  // nothing to see here
}
  
void EGS5AtmosphericCherenkovDetector::ausgab(int& iarg)
{
  if(iarg == 0)
    {
      extern EGS5_stack stack_;
      extern EGS5_epcont epcont_;
      const int ip = stack_.np-1;
      const int iq(stack_.iq[ip]);
      if(iq == 0)return;

      RadiatingTrackDetails rtd;
      rtd.e             = stack_.e[ip];
      rtd.x.set(stack_.x[ip], stack_.y[ip], stack_.z[ip]);
      rtd.t             = stack_.time[ip];
      rtd.u.set(stack_.u[ip], stack_.v[ip], stack_.w[ip]);
      rtd.ustep         = epcont_.ustep;
      rtd.n             = 1.0 + m_atm->nMinusOne(rtd.x.z());
      const double g2   = rtd.e*rtd.e/m_rmsq; // gamma^2
      const double b2   = 1.0 - 1.0/g2;       // beta^2
      rtd.sin2_thetac   = 1.0 - 1.0/(b2*rtd.n*rtd.n);
      if(rtd.sin2_thetac <= 0.0)goto skip_radiating_track;
      rtd.yield_density = 369.81020849958*rtd.sin2_thetac*rtd.ustep;
      rtd.cos_thetac    = std::sqrt(1.0 - rtd.sin2_thetac);
      rtd.sin_thetac    = std::sqrt(rtd.sin2_thetac);
      radiatingTrack(rtd);
    }
 skip_radiating_track: /* noop */ ;
}

void EGS5AtmosphericCherenkovDetector::
radiatingTrack(RadiatingTrackDetails& rtd)
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
//
// EGS5ACTArray
//
// ----------------------------------------------------------------------------

EGS5ACTArray::
EGS5ACTArray(Atmosphere& atm, unsigned nlayer, 
	     double emax, BField* bfield, double zbottom, double ztop,
	     unsigned nmedia)
: EGS5AtmosphericCherenkovDetector(atm, nlayer, emax, bfield, zbottom, ztop,
				   nmedia)
{
  // nothing to see here
}

EGS5ACTArray::~EGS5ACTArray()
{
  // nothing to see here
}
  
void EGS5ACTArray::radiatingTrack(RadiatingTrackDetails& rtd)
{
  for(std::vector<Scope>::const_iterator iscope = m_scopes.begin();
      iscope!=m_scopes.end(); iscope++)
    {
      Vec3D dx = iscope->x - rtd.x;
      double x2 = dx.norm2();
      double u  = dx*rtd.u;

      // Skip this telescope if particle is moving away from it
      if(u<0)continue; 

      double u2 = u*u;
      double v  = std::sqrt(x2-u2);
      double dmin = v*rtd.cos_thetac - u*rtd.sin_thetac;

      // Skip telescope if cone does not interesect with telescope sphere
      if(std::abs(dmin) > iscope->r)continue;

      double d2 = x2 - iscope->r*iscope->r;

      HitTelescopeDetails htd;
      htd.scope      = iscope;
      htd.iscope     = iscope-m_scopes.begin();
      htd.dxu        = u;
      htd.dxv        = v;
      htd.dmin       = dmin;
      htd.v          = (dx - rtd.u*u)/v;
      htd.cos_phimax = (std::sqrt(d2)-u*rtd.cos_thetac)/(v*rtd.sin_thetac);
      htd.phimax     = std::acos(std::max(htd.cos_phimax,-1.0));

      hitTelescope(htd, rtd);
    }
}

void EGS5ACTArray::hitTelescope(HitTelescopeDetails& htd, 
				RadiatingTrackDetails& rtd)
{
  // Default: do nothing
}

// ----------------------------------------------------------------------------
//
// EGS5SimpleIACTArray - Toy model of IACT without raytracing
//
// ----------------------------------------------------------------------------

void EGS5SimpleIACTArray::ImagingScope::
setZnAz(double zn, double az, double theta)
{
  RotationVec3D rvec(0,theta,0);
  rvec &= RotationVec3D(M_PI_2-zn,0,0);
  rvec &= RotationVec3D(0,0,-az);
  l = Vec3D(1,0,0); l.rotate(rvec);
  m = Vec3D(0,0,1); m.rotate(rvec);
  n = Vec3D(0,-1,0); n.rotate(rvec);
}

EGS5SimpleIACTArray::
EGS5SimpleIACTArray(Atmosphere& atm, unsigned nlayer, 
		    double emax, BField* bfield, double zbottom, double ztop,
		    unsigned nmedia)
: EGS5ACTArray(atm, nlayer, emax, bfield, zbottom, ztop, nmedia), m_iscopes()
{
  // nothing to see here
}

EGS5SimpleIACTArray::~EGS5SimpleIACTArray()
{
  // nothing to see here
}
  
void EGS5SimpleIACTArray::clearScopes() 
{
  m_iscopes.clear(); 
  EGS5ACTArray::clearScopes(); 
  m_images.clear();
}

void EGS5SimpleIACTArray::addScope(EGS5SimpleIACTArray::ImagingScope& s) 
{ 
  EGS5ACTArray::addScope(s); 
  m_iscopes.push_back(s); 
  m_images.push_back(Image(s.res, -0.5*s.fov, 0.5*s.fov));
}

void EGS5SimpleIACTArray::hitTelescope(HitTelescopeDetails& htd, 
				       RadiatingTrackDetails& rtd)
{
  const ImagingScope& scope(m_iscopes[htd.iscope]);
  Image& image(m_images[htd.iscope]);
  Vec3D ray0(rtd.cos_thetac*rtd.u+rtd.sin_thetac*htd.v);
  double observed_yield = rtd.yield_density*htd.phimax/M_PI;
  observed_yield *= scope.yield.yield(rtd.x.z(), -ray0.z());
  double half_arc_length = std::asin(rtd.sin_thetac)/M_PI*180.0*htd.phimax;
  int n = std::floor(5.0*half_arc_length/scope.res)+1;
  double dphi = half_arc_length/double(n);
  double ypp = observed_yield/double(2*n+1);

#if 1
  Vec3D w = rtd.u^htd.v;
  for(int i=-n;i<=n;i++)
    {
      double phi = double(i)*dphi;
      double cphi = std::cos(phi);
      double sphi = std::sin(phi);
      Vec3D ray(rtd.cos_thetac*rtd.u+rtd.sin_thetac*(htd.v*cphi+w*sphi));
      double l = (scope.l * ray)*180.0/M_PI;
      double m = (scope.m * ray)*180.0/M_PI;
      image.accumulate(l,m,ypp);
    }
#endif
}
