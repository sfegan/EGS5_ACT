// test.cpp - General test program for development
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: test.cpp 4954 2013-01-18 08:44:50Z sfegan $

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdint.h>
#include <map>

#include <VSOctaveH5Writer.hpp>

#include "EGS5LayeredDetector.hpp"
#include "EGS5AtmosphericDetector.hpp"
#include "DetectorEfficiency.hpp"
#include "Atmosphere.hpp"

#include "VSSimpleHist.hpp"
#include "RandomNumbers.hpp"

using namespace VERITAS;

// ****************************************************************************
// PREDEFINED SIMPLE RNG
// ****************************************************************************

class ConstantDeviateSimpleRNG: public SimpleRNG
{
public:
  ConstantDeviateSimpleRNG(double x): SimpleRNG(), m_x(x) { }
  virtual ~ConstantDeviateSimpleRNG();
  virtual double uniform();
private:
  double m_x;
};

ConstantDeviateSimpleRNG::~ConstantDeviateSimpleRNG()
{
  // nothing to see here
}

double ConstantDeviateSimpleRNG::uniform()
{
  return m_x;
}

class PredefinedDeviateSimpleRNG: public SimpleRNG
{
public:
  PredefinedDeviateSimpleRNG(SimpleRNG* rng);
  virtual ~PredefinedDeviateSimpleRNG();
  virtual double uniform();
  void addPD(uint64_t ideviate, double deviate);
  void resetSession();
  uint64_t sessionCount() { return m_session_count; }
  uint64_t totalCount() { return m_total_count+m_session_count; }
private:
  SimpleRNG*                           m_rng;
  uint64_t                             m_session_count;
  uint64_t                             m_total_count;
  std::map<uint64_t, double>           m_pd;
  std::map<uint64_t, double>::iterator m_ipd;
};

PredefinedDeviateSimpleRNG::~PredefinedDeviateSimpleRNG()
{
  // nothing to see here
}

double PredefinedDeviateSimpleRNG::uniform()
{
  double d;
  if(m_ipd!=m_pd.end() && m_session_count==m_ipd->first)d = (m_ipd++)->second;
  else d = m_rng->uniform();
  m_session_count++;
  return d;
}

PredefinedDeviateSimpleRNG::PredefinedDeviateSimpleRNG(SimpleRNG* rng):
  SimpleRNG(), m_rng(rng), m_session_count(), m_total_count(), m_pd(), 
  m_ipd(m_pd.begin())
{
  // nothing to see here  
}

void PredefinedDeviateSimpleRNG::addPD(uint64_t ideviate, double deviate)
{
  m_pd[ideviate]=deviate;
  resetSession();
}

void PredefinedDeviateSimpleRNG::resetSession() 
{
  m_total_count += m_session_count;
  m_session_count=0; 
  m_ipd = m_pd.begin();
}

// ****************************************************************************
// SLICE STAT DETECTOR
// ****************************************************************************

struct SliceStat
{
  SliceStat(double eres=100): 
    ntrack(), sumx(), sumy(), sumx2(), sumy2(), 
    e(), sumex(), sumey(), sumex2(), sumey2(), ex(eres,"ex"), ey(eres,"ey") { }
  unsigned ntrack;
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
private:
  unsigned               m_nevent;
  double                 m_slicesep;
  std::vector<SliceStat> m_slices;
  double                 m_eres;
};

SliceStatDetector::
SliceStatDetector(const std::vector<Media>& media, const std::vector<Layer>& layers, 
	     double zbot, double emax, double slicesep, double eres):
  EGS5LayeredDetector(media, layers, zbot, emax), 
  m_nevent(), m_slicesep(slicesep), m_slices(0,eres), m_eres(eres)
{
  // nothing to see here
}

SliceStatDetector::~SliceStatDetector()
{
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
	std::cout 
	  << islice << ' ' 
	  << slice.ntrack << ' ' << meanx << ' ' << meany << ' ' 
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
// TRACK WRITING DETECTOR
// ****************************************************************************

class TrackWritingDetector: public EGS5LayeredDetector
{
public:
  TrackWritingDetector(const std::vector<Media>& media, 
		       const std::vector<Layer>& layers, 
		       double zbot, double emax);
  virtual ~TrackWritingDetector();
  virtual void showerCompleted();
  virtual void ausgab(int& iarg);  
private:
  unsigned               m_nevent;
};

TrackWritingDetector::
TrackWritingDetector(const std::vector<Media>& media, 
		     const std::vector<Layer>& layers, 
		     double zbot, double emax):
  EGS5LayeredDetector(media, layers, zbot, emax), m_nevent()
{
  // nothing to see here
}

TrackWritingDetector::~TrackWritingDetector()
{
  // nothing to see here
}

void TrackWritingDetector::showerCompleted()
{
  m_nevent++;
}

void TrackWritingDetector::ausgab(int& iarg)
{
  extern EGS5_stack stack_;
  extern EGS5_epcont epcont_;
  int ip = stack_.np-1;

  if(iarg == 0)
    {
      std::cout 
	<< m_nevent << ' ' << stack_.e[ip] << ' '  << stack_.iq[ip] << ' ' 
	<< stack_.time[ip] << ' '
	<< stack_.x[ip] << ' ' << stack_.y[ip] << ' ' << stack_.z[ip] << ' '
	<< stack_.x[ip]+stack_.u[ip]*epcont_.ustep << ' '
	<< stack_.y[ip]+stack_.v[ip]*epcont_.ustep << ' '
	<< stack_.z[ip]+stack_.w[ip]*epcont_.ustep << '\n';
    }
}

// ****************************************************************************
// MAIN
// ****************************************************************************

typedef EGS5LayeredDetector::Media Media;
typedef EGS5LayeredDetector::Layer Layer;

int main()
{
  TelescopeEfficiency eff;
  eff.scaleEffFromFile("Parameters/corsika_mirreff.dat");
  eff.scaleEffFromFile("Parameters/corsika_quanteff.dat");
  AtmosphericAbsorption atmabs("Parameters/corsika_atmabs.dat");
#if 0
  for(double d=0; d<6000000; d+=10000)
    {
      yield_t y = yield(d);
      std::cout << d/100000 << ' '
		<< y.first << ' ' << y.second << ' ' << y.third << '\n';
    }
  std::exit(1);
#endif
  try
    {
      SimpleRNG* base_rng = EGS5RanluxSimpleRNG::instance();

      //RandomNumbers rng_base(RandomNumbers::defaultFilename());
      //SimpleRNGAdapter rng(&rng_base);
      //base_rng = &rng;

      PredefinedDeviateSimpleRNG pd_rng(base_rng);
      pd_rng.addPD(0,0.3);

      ConstantDeviateSimpleRNG c_rng(0.5);

      LayeredAtmosphere atm("Parameters/atmprof6.dat");
      std::vector<Media> media;
      std::vector<Layer> layers;
      double zbottom = 0;
      double ztop = HUGE_VAL;
      ConstantBField const_bfield(  -2831.9 * 1e-5,  
				    11719.0 * 1e-5,
				   -25688.5 * 1e-5);
      BField* bfield = &const_bfield;
      EGS5AtmosphericDetector::
	makeMediaAndLayers(media, layers, atm, 100, bfield, zbottom, ztop);

      //EGS5LayeredDetector det(media, layers, zbottom, 1000000);

      EGS5SimpleIACTArray det(atm, 100, 1000000, bfield, zbottom, ztop, 1);
      double z = 150000.0;
      EGS5SimpleIACTArray::ImagingScope s; 
      s.res = 0.02;
      s.fov = 5.0;
      s.x.set( 6000, 6000, z); s.r=600;  det.addScope(s);
      s.x.set(-6000, 6000, z); s.r=600;  det.addScope(s);
      s.x.set( 6000,-6000, z); s.r=600;  det.addScope(s);
      s.x.set(-6000,-6000, z); s.r=600;  det.addScope(s);
      s.x.set(    0,    0, z); s.r=1400; det.addScope(s);

      double theta = 0/180.0*M_PI;
      double ctheta = std::cos(theta);
      double stheta = std::sin(theta);

      for(std::vector<EGS5SimpleIACTArray::ImagingScope>::iterator
	    iscope = det.scopes().begin(); iscope != det.scopes().end();
	  iscope++)
	{
	  iscope->setZnAz(theta, 0);
	  iscope->yield = atmabs.integrateYield(iscope->x.z(), ctheta, eff);
	}

      EGS5System* egs5 = EGS5System::instance(&det, base_rng);
      //egs5->setRNG(35, &c_rng); // Number of MFP for hard step in e+/e-
      //egs5->setRNG(86, &pd_rng); // Number of MFP for photon      
      egs5->initializeEGS5();
      
      for(unsigned i=0;i<100;i++)
	{
	  pd_rng.resetSession();
	  egs5->shower(100000, 6000, 6000, layers.back().zt,
		       stheta, 0.0, -ctheta,
		       layers.size()+1, 0);
	  std::cout << i+1 << '\n';
	}

      unsigned nscope = det.images().size();
      VSOctaveH5Writer writer("doodles/image.h5");
      VSOctaveH5WriterCellVector* c = writer.writeCellVector("images",nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  VSOctaveH5WriterStruct* s = c->writeStruct(iscope);
	  det.images()[iscope].save(s);
	  delete s;
	}
      delete c;
    }
  catch(std::string& s)
    {
      std::cerr << s << '\n';
    }
  catch(const char* s)
    {
      std::cerr << s << '\n';
    }
  std::cerr << "Fin...\n";
}
