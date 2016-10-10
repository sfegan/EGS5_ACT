// SimpleInstrumentedDetectors.cpp - Some simple detectors to measure
// shower profiles and make track diagrams
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5LayeredDetector.cpp 4954 2013-01-18 08:44:50Z sfegan $

#include <fstream>
#include <iostream>

#include "SimpleInstrumentedDetectors.hpp"

// ****************************************************************************
// TRACK COUNTING DETECTOR
// ****************************************************************************

TrackCountingDetector::TrackCountingDetector():
  EGS5UserInterface(), m_ph_ntrack(), m_ep_ntrack()
{
  // nothing to see here
}

TrackCountingDetector::~TrackCountingDetector()
{
  // nothing to see here
}

double TrackCountingDetector::getEMax()
{
  return 0;
}

void TrackCountingDetector::getMedia(std::vector<Media>& media)
{
  // nothing to see here
}

void TrackCountingDetector::getRegions(std::vector<Region>& regions)
{
  // nothing to see here
}

void TrackCountingDetector::howfar()
{
  // nothing to see here
}

void TrackCountingDetector::ausgab(int& iarg)
{
  if(iarg == 0)
    {
      extern EGS5_stack stack_;
      extern EGS5_epcont epcont_;
      int ip = stack_.np-1;
      if(stack_.iq[ip] == 0)
	m_ph_ntrack++;
      else
	m_ep_ntrack++;
    }
}

// ****************************************************************************
// TRACK WRITING DETECTOR
// ****************************************************************************

TrackWritingDetector::
TrackWritingDetector(const std::vector<Media>& media, 
		     const std::vector<Layer>& layers, 
		     double zbot, double emax, double ecut):
  EGS5LayeredDetector(media, layers, zbot, emax), m_tracks(), m_ecut(ecut)
{
  // nothing to see here
}

TrackWritingDetector::~TrackWritingDetector()
{
  // nothing to see here
}

void TrackWritingDetector::
showerStarting(double e,
	       double x, double y, double z,
	       double u, double v, double w,
	       int iregion, int q, double weight)
{
  clear();
}

void TrackWritingDetector::ausgab(int& iarg)
{
  if(iarg == 0)
    {
      extern EGS5_stack stack_;
      extern EGS5_epcont epcont_;
      int ip = stack_.np-1;
      if(stack_.e[ip] > m_ecut)
	{
	  Track t;
	  t.e   = stack_.e[ip];
	  t.iq  = stack_.iq[ip];
	  t.t0  = stack_.time[ip];
	  t.x0  = stack_.x[ip];
	  t.y0  = stack_.y[ip];
	  t.z0  = stack_.z[ip];
	  t.x1  = stack_.x[ip]+stack_.u[ip]*epcont_.ustep;
	  t.y1  = stack_.y[ip]+stack_.v[ip]*epcont_.ustep;
	  t.z1  = stack_.z[ip]+stack_.w[ip]*epcont_.ustep;
	  m_tracks.push_back(t);
	}
    }
}

// ****************************************************************************
// SliceStat Detector
// ****************************************************************************

SliceStatDetector::
SliceStatDetector(const std::vector<Media>& media, 
		  const std::vector<Layer>& layers, 
		  double zbot, double emax, double slicesep, double eres):
  EGS5LayeredDetector(media, layers, zbot, emax), 
  m_nevent(), m_slicesep(slicesep), m_slices(0,eres), m_eres(eres)
{
  // nothing to see here
}

SliceStatDetector::~SliceStatDetector()
{

}

void SliceStatDetector::dump(const std::string filename)
{
  std::ostream *str = &std::cout;
  std::ofstream fstr;

  if(filename != std::string())
    {
      fstr.open(filename.c_str());
      str = &fstr;
    }

  const double efrac = 0.8;
  for(unsigned islice=0;islice<m_slices.size();islice++)
    if(m_slices[islice].ntrack)
      {
	SliceStat& slice(m_slices[islice]);
	double ntrack = double(slice.ntrack);
	double meanx = slice.sumx/ntrack;
	double meany = slice.sumy/ntrack;
	double varx = slice.sumx2/ntrack - meanx*meanx;
	double vary = slice.sumy2/ntrack - meany*meany;
	double e = slice.e;
	double meanex = slice.sumex/e;
	double meaney = slice.sumey/e;
	double varex = slice.sumex2/e - meanex*meanex;
	double varey = slice.sumey2/e - meaney*meaney;

	VSSimpleHist<double,double> iex = slice.ex.makeNormIntegralHist();
	double exw = iex.integralCountToVal(efrac);
	VSSimpleHist<double,double> iey = slice.ey.makeNormIntegralHist();
	double eyw = iey.integralCountToVal(efrac);
	(*str)
	  << islice << ' ' 
	  << slice.ntrack << ' ' 
	  << slice.ncharged << ' ' << slice.nphoton << ' ' 
	  << meanx << ' ' << meany << ' ' 
	  << std::sqrt(varx) << ' ' << std::sqrt(vary) << ' '
	  << e << ' ' << meanex << ' ' << meany << ' ' 
	  << std::sqrt(varex) << ' ' << std::sqrt(varey) << ' ' 
	  << exw << ' ' << eyw << '\n';
      }
}

void SliceStatDetector::
showerStarting(double e, double x, double y, double z,
	       double u, double v, double w, int iregion, int q, double weight)
{
#if 0
  for(unsigned islice=0;islice<m_slices.size();islice++)
    m_slices[islice]=SliceStat();
#endif
  m_nevent++;
}

void SliceStatDetector::showerCompleted()
{
  
}

void SliceStatDetector::ausgab(int& iarg)
{
  EGS5LayeredDetector::ausgab(iarg);
  
  extern EGS5_stack stack_;
  extern EGS5_epcont epcont_;
  int ip = stack_.np-1;

  if(iarg == 0)
    {
      double xnow = stack_.x[ip];
      double ynow = stack_.y[ip];
      double znow = stack_.z[ip];
      int islicenow = std::floor(znow / m_slicesep);
      double znxt = znow + stack_.w[ip]*epcont_.ustep;
      int islicenxt = std::floor(znxt / m_slicesep);
      islicenow = std::min(islicenow, 1500);
      islicenow = std::max(islicenow, 0);
      if(islicenow >= m_slices.size())
	m_slices.resize(islicenow+1, SliceStat(m_eres));
      while((islicenow>=0)&&(islicenow>islicenxt))
	{
	  SliceStat& slice(m_slices[islicenow]);
	  double zslice = double(islicenow)*m_slicesep;
	  double xslice = xnow + stack_.u[ip]/stack_.w[ip]*(zslice-znow);
	  double yslice = ynow + stack_.v[ip]/stack_.w[ip]*(zslice-znow);
	  double e = stack_.e[ip];
	  slice.ntrack++;
	  if(stack_.iq[ip] == 0)
	    slice.nphoton++;
	  else
	    slice.ncharged++;
	  slice.sumx   += xslice;
	  slice.sumy   += yslice;
	  slice.sumx2  += xslice*xslice;
	  slice.sumy2  += yslice*yslice;
	  slice.e      += e;
	  slice.sumex  += e*xslice;
	  slice.sumey  += e*yslice;
	  slice.sumex2 += e*xslice*xslice;
	  slice.sumey2 += e*yslice*yslice;
	  slice.ex.accumulate(std::min(std::fabs(xslice),100000.0),e);
	  slice.ey.accumulate(std::min(std::fabs(yslice),100000.0),e);
	  islicenow--;
	}
    }
}

// ****************************************************************************
// Total Cherenkov Yield on Ground
// ****************************************************************************

TotalCherenkovYieldDetector::
TotalCherenkovYieldDetector(Atmosphere& atm, unsigned nlayer, double emax, 
			    ACTIntegratedLightYield* ground_yield,
			    BField* bfield, 
			    double zbottom, double ztop,
			    unsigned nmedia):
  EGS5AtmosphericCherenkovDetector(atm, nlayer, emax, bfield, zbottom, ztop,
				   nmedia),
  m_ground_yield(ground_yield), m_yield(0)
{

}

void TotalCherenkovYieldDetector::radiatingTrack(RadiatingTrackDetails& rtd)
{
  double observed_yield = rtd.yield_density;
  observed_yield *= m_ground_yield->yield(rtd.x.z(), -rtd.u.z());
  m_yield += observed_yield;
}

// ****************************************************************************
// PRUNING DETECTOR
// ****************************************************************************

PruningDetector::PruningDetector(double ecut):
  EGS5UserInterface(), m_particles(), m_ecut(ecut)
{
  // nothing to see here
}

PruningDetector::~PruningDetector()

{
  // nothing to see here
}

void PruningDetector::showerStarting(double e,
				     double x, double y, double z,
				     double u, double v, double w,
				     int iregion, int q, double weight)
{
  clear();
}

double PruningDetector::getEMax()
{
  return 0;
}

void PruningDetector::getMedia(std::vector<Media>& media)
{
  // nothing to see here
}

void PruningDetector::getRegions(std::vector<Region>& regions)
{
  // nothing to see here
}

void PruningDetector::ausgab(int& iarg)
{
  // nothing to see here
}

void PruningDetector::howfar()
{
  extern EGS5_stack stack_;
  extern EGS5_epcont epcont_;

  int ip = stack_.np-1;
  double e = stack_.e[ip];

  if(e <= m_ecut)
    {
      Particle p;
      p.e  = e;
      p.iq = stack_.iq[ip];
      p.x  = stack_.x[ip];
      p.y  = stack_.y[ip];
      p.z  = stack_.z[ip];
      p.u  = stack_.u[ip];
      p.v  = stack_.v[ip];
      p.w  = stack_.w[ip];
      p.t  = stack_.time[ip];
      p.ir = stack_.ir[ip]-1;
      m_particles.push_back(p);
      epcont_.idisc = 1;
      return;
    }
}
