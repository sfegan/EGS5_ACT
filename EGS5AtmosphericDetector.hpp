// EGS5AtmosphericDetector.hpp - Detector class for EGS5 system that
// - implements a layered atmospheric detector
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2012
// $Id: EGS5AtmosphericDetector.hpp 5688 2013-09-30 13:58:43Z sfegan $

#ifndef EGS5ATMOSPHERICDETECTOR_HPP
#define EGS5ATMOSPHERICDETECTOR_HPP

#include <stdint.h>
#include <cmath>

#include <VSAAlgebra.hpp>
#include <VSSimple2DHist.hpp>

typedef VERITAS::VSAAlgebra::Vec3D               Vec3D;
typedef VERITAS::VSAAlgebra::RotationVec3D       RotationVec3D;

#include "EGS5LayeredDetector.hpp"
#include "Atmosphere.hpp"
#include "BField.hpp"
#include "DetectorEfficiency.hpp"

// ****************************************************************************
// EGS5AtmosphericDetector - convenience class to construct layered detector
// from atmosphere and magnetic field
// ****************************************************************************

class EGS5AtmosphericDetector: public EGS5LayeredDetector
{
public:
  EGS5AtmosphericDetector(Atmosphere& atm, unsigned nlayer, double emax, 
			  BField* bfield = 0, 
			  double zbottom = 0, double ztop = HUGE_VAL,
			  unsigned nmedia = 1);
  virtual ~EGS5AtmosphericDetector();
  
  virtual void writePEGS5InputFile(const std::string& filename);

  static void makeMediaAndLayers(std::vector<Media>& media,
				 std::vector<Layer>& layers,
				 Atmosphere& atm, unsigned nlayer,
				 BField* bfield = 0,
				 double zbottom = 0, double ztop = HUGE_VAL,
				 unsigned nmedia = 1);

protected:
  void initialize(Atmosphere* atm, unsigned nlayer, double emax,
		  BField* bfield, double zbottom, double ztop,
		  unsigned nmedia);

  Atmosphere*         m_atm;
  BField*             m_bfield;
  std::vector<std::pair<Media,double> > m_air_media;
};

// ****************************************************************************
// EGS5AtmosphericCherenkovDetector - Test whether charged particles are
// above the Cherenkov threshold while being transported.
// ****************************************************************************

class EGS5AtmosphericCherenkovDetector: public EGS5AtmosphericDetector
{
public:
  EGS5AtmosphericCherenkovDetector(Atmosphere& atm, unsigned nlayer, 
				   double emax, BField* bfield = 0, 
				   double zbottom = 0, double ztop = HUGE_VAL,
				   unsigned nmedia = 1);
  virtual ~EGS5AtmosphericCherenkovDetector();
  
  virtual void ausgab(int& iarg);  

  class RadiatingTrackDetails
  {
  public:
    double sin2_thetac;    // Cosine^2 of emission angle   [1]
    double sin_thetac;     // Sine of emssion angle        [1]
    double cos_thetac;     // Cosine of emission angle     [1]
    double yield_density;  // Cherenkov photon density     [ph/eV]
    double n;              // Refracive index              [1]
    double e;              // Particle energy              [MeV]
    Vec3D x;               // Position of track            [cm]
    double t;              // Time at start of track       [ns]
    Vec3D u;               // Direction of  propogation    [1]
    double ustep;          // Track length                 [cm]
  };

  virtual void radiatingTrack(RadiatingTrackDetails& rtd);
};

// ****************************************************************************
// EGS5ACTArray - Propagate photon yield to telescopes on ground
// ****************************************************************************

class EGS5ACTArrayScope
{
public:
  Vec3D x;
  double r;
};

class EGS5ACTArray: public EGS5AtmosphericCherenkovDetector
{
public:
  typedef EGS5ACTArrayScope Scope;

  class HitTelescopeDetails
  {
  public:
    std::vector<Scope>::const_iterator scope;
    unsigned iscope;       //
    Vec3D v;               //
    double dxu;            //
    double dxv;            //
    double dmin;           //
    double cos_phimax;     //
    double phimax;         //
  };

  EGS5ACTArray(Atmosphere& atm, unsigned nlayer, 
	       double emax, BField* bfield = 0, 
	       double zbottom = 0, double ztop = HUGE_VAL,
	       unsigned nmedia = 1);
  virtual ~EGS5ACTArray();
  
  void clearScopes() { m_scopes.clear(); }
  std::vector<Scope>& scopes() { return m_scopes; }
  void addScope(Scope& s) { m_scopes.push_back(s); }

  virtual void radiatingTrack(RadiatingTrackDetails& rtd);
  virtual void hitTelescope(HitTelescopeDetails& htd, 
			    RadiatingTrackDetails& rtd);

protected:
  std::vector<Scope> m_scopes;
};

// ****************************************************************************
// EGS5SimpleIACTArray - Toy model of IACT without raytracing
// ****************************************************************************

class EGS5ACTArrayImagingScope: public EGS5ACTArrayScope
{
public:
  double w0; // Focus parameter: 1/image distance telescope is focused at
  Vec3D l; // Vector along "x-axis" of imaging plane
  Vec3D m; // Vector along "y-axis" of imaging plane
  Vec3D n; // Vector along telescope pointing direction
  double res;
  double fov;
  ACTIntegratedLightYield yield;
  void setZnAz(double zn, double az, double theta=0);
};

typedef VERITAS::VSSimple2DHist<double, double> EGS5ACTImage;

class EGS5SimpleIACTArray: public EGS5ACTArray
{
public:
  typedef EGS5ACTArrayImagingScope ImagingScope;
  typedef EGS5ACTImage             Image;

  EGS5SimpleIACTArray(Atmosphere& atm, unsigned nlayer, 
		      double emax, BField* bfield = 0, 
		      double zbottom = 0, double ztop = HUGE_VAL,
		      unsigned nmedia = 1);
  virtual ~EGS5SimpleIACTArray();
  
  void clearScopes();
  std::vector<ImagingScope>& scopes() { return m_iscopes; }
  void addScope(ImagingScope& s); 

  const Image& image(unsigned iimage) const { return m_images[iimage]; }
  const std::vector<Image>& images() const { return m_images; }

  Image& image(unsigned iimage) { return m_images[iimage]; }
  std::vector<Image>& images() { return m_images; }
  void clearImages() { for(unsigned iimage=0;iimage<m_images.size();iimage++)m_images[iimage].clear(); }

  virtual void hitTelescope(HitTelescopeDetails& htd, 
			    RadiatingTrackDetails& rtd);

protected:
  std::vector<ImagingScope>   m_iscopes;
  std::vector<Image> m_images;
};

#endif // define EGS5ATMOSPHERICDETECTOR_HPP
