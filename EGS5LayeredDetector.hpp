// EGS5LayeredDetector.hpp - Detector class for EGS5 system that implements
// - layered detector geometry and magnetic field (as in ucbend.f).
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5LayeredDetector.hpp 5422 2013-06-26 14:01:03Z sfegan $

#ifndef EGS5LAYEREDDETECTOR_HPP
#define EGS5LAYEREDDETECTOR_HPP

#include <stdint.h>
#include <cmath>

#include "EGS5System.hpp"

// Class describing a "sandwich" detector made up of infinite parallel
// layers stacked in the z-direction, optionally with a magnetic field
// in each region. The detector sits in a vacuum for z less than some
// zbot and for z greater than

struct EGS5Layer
{
  EGS5Layer(double _zt=0):
    imedia(), rho(), zt(_zt), bfield(), bx(), by(), bz() { }
  int imedia;    // Index into media array - (-1 is vacuum)
  double rho;    // Override default media density
  double zt;     // "Top" of region on +z direction
  bool   bfield; // Magnetic field enabled in region;
  double bx;     // X-component of b-field
  double by;     // Y-component of b-field
  double bz;     // Z-component of b-field
};

class EGS5LayeredDetector: public EGS5UserInterface
{
public:
  typedef EGS5Layer Layer;

  EGS5LayeredDetector(const std::vector<Media>& media, 
		      const std::vector<Layer>& layers, 
		      double zbot, double emax);
  virtual ~EGS5LayeredDetector();

  inline unsigned getLayerNumber(double z) const;

  // **************************************************************************
  // REQUIRED FUNCTIONS - must be supplied by derived class
  // **************************************************************************

  // Return maximum energy
  virtual double getEMax();

  // Return list of media names
  virtual void getMedia(std::vector<Media>& media);

  // Return list of regions in problem
  virtual void getRegions(std::vector<Region>& regions);

  // Primary EGS5 interface
  virtual void ausgab(int& iarg);
  virtual void howfar(void);

  // Return vector or layers
  void getLayers(std::vector<Layer>& layers);

protected:
  EGS5LayeredDetector();
  void initialize(const std::vector<Media>& media, 
		  const std::vector<Layer>& layers, 
		  double zbot, double emax);

  struct InternalLayer: public Layer
  {
    InternalLayer(const Layer& _l = Layer(), double _zb = 0): 
      Layer(_l), zb(_zb), 
      btotal(std::sqrt(_l.bx*_l.bx + _l.by*_l.by + _l.bz*_l.bz)),
      ubx(btotal>0?_l.bx/btotal:0), uby(btotal>0?_l.by/btotal:0),
      ubz(btotal>0?_l.bz/btotal:0) { }
    double zb;
    double btotal;
    double ubx;
    double uby;
    double ubz;
    bool operator< (const InternalLayer& o) const { return zt<o.zt;  }
  };

  std::vector<Media>         m_media;
  std::vector<InternalLayer> m_layers;
  double                     m_zbot;
  double                     m_emax;
  double                     m_rmsq;
  uint64_t                   m_n_howfar_entry;
  uint64_t                   m_n_howfar_change_region;
};

unsigned EGS5LayeredDetector::getLayerNumber(double z) const
{
  if(z<m_zbot)return 0;
  std::vector<InternalLayer>::const_iterator ilayer =
    std::lower_bound(m_layers.begin(), m_layers.end(), InternalLayer(z));
  return (ilayer-m_layers.begin())+1;
}

#endif // not defined EGS5LAYEREDDETECTOR_HPP
