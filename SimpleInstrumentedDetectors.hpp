// SimpleInstrumentedDetectors.hpp - Some simple detectors to measure
// shower profiles and make track diagrams
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5LayeredDetector.cpp 4954 2013-01-18 08:44:50Z sfegan $

#ifndef SIMPLEINSTRUMENTEDDETECTORS_HPP
#define SIMPLEINSTRUMENTEDDETECTORS_HPP

#include "VSSimpleHist.hpp"

#include "EGS5LayeredDetector.hpp"
#include "EGS5AtmosphericDetector.hpp"

using namespace VERITAS;

// ****************************************************************************
// TRACK COUNTING DETECTOR
// ****************************************************************************

class TrackCountingDetector: public EGS5UserInterface
{
public:
  TrackCountingDetector();
  virtual ~TrackCountingDetector();

  virtual double getEMax();
  virtual void getMedia(std::vector<Media>& media);
  virtual void getRegions(std::vector<Region>& regions);
  virtual void ausgab(int& iarg);  
  virtual void howfar();  

  unsigned nTracks() const { return m_ph_ntrack+m_ep_ntrack; }
  unsigned nPhotonTracks() const { return m_ph_ntrack; }
  unsigned nChargedTracks() const { return m_ep_ntrack; }

  void reset() { m_ph_ntrack = m_ep_ntrack = 0; }
private:
  unsigned           m_ph_ntrack;
  unsigned           m_ep_ntrack;
};

// ****************************************************************************
// TRACK WRITING DETECTOR
// ****************************************************************************

struct Track
{
  double e;
  int    iq;
  double t0;
  double x0;
  double y0;
  double z0;
  double x1;
  double y1;
  double z1;
};

class TrackWritingDetector: public EGS5LayeredDetector
{
public:
  TrackWritingDetector(const std::vector<Media>& media, 
		       const std::vector<Layer>& layers, 
		       double zbot, double emax,
		       double ecut = 0);
  virtual ~TrackWritingDetector();
  virtual void showerStarting(double e,
			      double x, double y, double z,
			      double u, double v, double w,
			      int iregion, int q, double weight);
  virtual void ausgab(int& iarg);  

  unsigned nTracks() const { return m_tracks.size(); }
  const Track& track(unsigned itrack) const { return m_tracks[itrack]; }
  const std::vector<Track>& tracks() const { return m_tracks; }
  void clear() { m_tracks.clear(); }
private:
  std::vector<Track> m_tracks;
  double             m_ecut;
};

// ****************************************************************************
// SliceStat Detector
// ****************************************************************************

struct SliceStat
{
  SliceStat(double eres=100): 
    ntrack(), ncharged(), nphoton(), sumx(), sumy(), sumx2(), sumy2(), 
    e(), sumex(), sumey(), sumex2(), sumey2(), ex(eres,"ex"), ey(eres,"ey") { }
  unsigned ntrack;
  unsigned ncharged;
  unsigned nphoton;
  double sumx;
  double sumy;
  double sumx2;
  double sumy2;
  double e;
  double sumex;
  double sumey;
  double sumex2;
  double sumey2;
  VSSimpleHist<double,double> ex;
  VSSimpleHist<double,double> ey;  
};

class SliceStatDetector: public EGS5LayeredDetector
{
public:
  SliceStatDetector(const std::vector<Media>& media, 
		    const std::vector<Layer>& layers, 
		    double zbot, double emax, double slicesep=10000.0, 
		    double eres=100);
  virtual ~SliceStatDetector();
  virtual void showerStarting(double e,
			      double x, double y, double z,
			      double u, double v, double w,
			      int iregion, int q, double weight);
  virtual void showerCompleted();
  virtual void ausgab(int& iarg);  
  void dump(const std::string filename = std::string());
private:
  unsigned               m_nevent;
  double                 m_slicesep;
  std::vector<SliceStat> m_slices;
  double                 m_eres;
};

// ****************************************************************************
// Total Cherenkov Yield on Ground
// ****************************************************************************

class TotalCherenkovYieldDetector: public EGS5AtmosphericCherenkovDetector
{
public:
  TotalCherenkovYieldDetector(Atmosphere& atm, unsigned nlayer, double emax, 
			      ACTIntegratedLightYield* ground_yield,
			      BField* bfield = 0, 
			      double zbottom = 0, double ztop = HUGE_VAL,
			      unsigned nmedia = 1);
  virtual void radiatingTrack(RadiatingTrackDetails& rtd);
  void reset() { m_yield = 0; }
  double yield() const { return m_yield; }

private:
  ACTIntegratedLightYield* m_ground_yield;
  double m_yield;
};

// ****************************************************************************
// Pruning Detector
// ****************************************************************************

struct Particle
{
  double e;
  int    iq;
  double x;
  double y;
  double z;
  double u;
  double v;
  double w;
  double t;
  int    ir;  // EGS5 region number - 1 (i.e. starts from 0)
};

class PruningDetector: public EGS5UserInterface
{
public:
  PruningDetector(double ecut);
  virtual ~PruningDetector();
  virtual void showerStarting(double e,
			      double x, double y, double z,
			      double u, double v, double w,
			      int iregion, int q, double weight);

  virtual double getEMax();
  virtual void getMedia(std::vector<Media>& media);
  virtual void getRegions(std::vector<Region>& regions);
  virtual void ausgab(int& iarg);  
  virtual void howfar();  
  
  unsigned nParticles() const { return m_particles.size(); }
  const Particle& particle(unsigned ip) const { return m_particles[ip]; }
  const std::vector<Particle>& particles() const { return m_particles; }
  void clear() { m_particles.clear(); }
private:
  std::vector<Particle> m_particles;
  double                m_ecut;
};

// ****************************************************************************
// Pruning Detector
// ****************************************************************************


#endif // defined SIMPLEINSTRUMENTEDDETECTORS_HPP
