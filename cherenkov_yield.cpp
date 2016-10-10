// cherenkov_yield.cpp - Calculate the mean photon yield distribution
// (and profile with heigth) from a primary of predefined energy and
// start height.
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012 
// $Id: cherenkov_yield.cpp 4591 2012-10-08 14:00:50Z sfegan $

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

#include "EGS5LayeredDetector.hpp"
#include "EGS5AtmosphericDetector.hpp"
#include "Atmosphere.hpp"

#include "VSSimpleHist.hpp"
#include "VSOptions.hpp"
#include "RandomNumbers.hpp"

#define MELEC 0.51099891

using namespace VERITAS;

class ACDYield: public EGS5AtmosphericCherenkovDetector
{
public:
  ACDYield(Atmosphere& atm, unsigned nlayer, double emax, BField* bfield, 
	   double yield_res = 1.0, double zres = 10000.0,
	   double zbottom = 0, double ztop = HUGE_VAL, unsigned nmedia = 1);
  virtual ~ACDYield();
  virtual void radiatingTrack(double sin2_thetac, double n, double e, 
			      double x, double y, double z,
			      double ux, double uy, double uz,
			      double ustep);
  virtual void showerCompleted();

  void yieldStat(double& mean, double& rms);

  void reset();
  void writeResults(VSOctaveH5WriterStruct* s);

private:
  typedef VSSimpleHist<double,double> ProfileHist;
  typedef VSSimpleHist<double,unsigned> YieldHist;
  unsigned    m_nshower;
  double      m_total_yield;
  double      m_total_yield_sq;
  YieldHist   m_shower_yield_hist;
  ProfileHist m_shower_yield_profile;
  ProfileHist m_total_yield_profile;
  ProfileHist m_total_yield_sq_profile;
};

ACDYield::
ACDYield(Atmosphere& atm, unsigned nlayer, double emax, BField* bfield, 
	 double yield_res, double zres, double zbottom, double ztop, 
	 unsigned nmedia):
  EGS5AtmosphericCherenkovDetector(atm, nlayer, emax, bfield, zbottom, ztop,
				   nmedia),
  m_nshower(), m_total_yield(), m_total_yield_sq(),
  m_shower_yield_hist(yield_res), m_shower_yield_profile(zres), 
  m_total_yield_profile(zres), m_total_yield_sq_profile(zres)
{
  // nothing to see here
}

void ACDYield::reset()
{
  m_nshower        = 0;
  m_total_yield    = 0;
  m_total_yield_sq = 0;
  m_shower_yield_hist.clear();
  m_shower_yield_profile.clear();
  m_total_yield_profile.clear();
  m_total_yield_sq_profile.clear();
}

ACDYield::~ACDYield()
{
  // nothing to see here
}

void ACDYield::radiatingTrack(double sin2_thetac, double n, double e, 
			      double x, double y, double z,
			      double ux, double uy, double uz,
			      double ustep)
{
  double yield_per_ev = 369.81020849958*sin2_thetac*ustep;
#if 0
  double ethresh = std::sqrt(m_rmsq / (1.0 - 1.0/(n*n)));
  double thetac = std::asin(std::sqrt(sin2_thetac))*180.0/M_PI;
  std::cout 
    << e << ' ' << z/100000 << ' '
    << n << ' ' << thetac << ' ' << ethresh << ' ' << yield_per_ev << '\n';
#endif
  m_shower_yield_profile.accumulate(z, yield_per_ev);
}

void ACDYield::showerCompleted()
{
  m_nshower++;
  double shower_yield = m_shower_yield_profile.sum();
  m_total_yield += shower_yield;
  m_total_yield_sq += shower_yield*shower_yield;
  m_shower_yield_hist.accumulate(shower_yield);
  m_total_yield_profile.injest(m_shower_yield_profile, ProfileHist::CountOp());
  m_total_yield_sq_profile.injest(m_shower_yield_profile, ProfileHist::CountSqOp());
  m_shower_yield_profile.clear();
}

void ACDYield::yieldStat(double& mean, double& rms)
{
  mean = m_total_yield/double(m_nshower);
  rms = std::sqrt(m_total_yield_sq/double(m_nshower)-mean*mean);
}

void ACDYield::
writeResults(VSOctaveH5WriterStruct* s)
{
  VSOctaveH5WriterStruct* s2;
  VSOctaveH5WriterStruct* s3;
  s->writeScalar("nshower", m_nshower);

  double mean_yield, rms_yield;
  yieldStat(mean_yield, rms_yield);
  s2 = s->writeStruct("total_yield");
  s2->writeScalar("mean", mean_yield);
  s2->writeScalar("std", rms_yield);
  delete s2;

  s2 = s->writeStruct("profile");
  s3 = s2->writeStruct("mean");
  ProfileHist prof_mean(m_total_yield_profile);
  prof_mean /= m_nshower;
  prof_mean.save(s3);
  delete s3;
  ProfileHist prof_std(m_total_yield_sq_profile);
  prof_std  /= m_nshower;
  prof_mean *= prof_mean;
  prof_std  -= prof_mean;
  prof_std.transform(ProfileHist::CountSqrtOp());
  s3 = s2->writeStruct("std");
  prof_std.save(s3);
  delete s3;
  delete s2;

  s2 = s->writeStruct("shower_yield");
  m_shower_yield_hist.save(s2);
  delete s2;

  std::cout << mean_yield << ' ' << rms_yield << std::endl;
}

// ****************************************************************************
// MAIN
// ****************************************************************************

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] command_argument"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

typedef EGS5LayeredDetector::Media Media;
typedef EGS5LayeredDetector::Layer Layer;

double cherenkovThreshold(Atmosphere& atm, double h)
{
  double n = atm.nMinusOne(h)+1.0;
  double eth = MELEC/std::sqrt(1.0-1.0/(n*n));
  return eth;
}

void runSim(double h, double e, double ux, double uy, double uz, int q,
	    Atmosphere& atm, EGS5System* egs5,
	    ACDYield& det,  unsigned nsim, VSOctaveH5WriterStruct* s)
{
  double eth = cherenkovThreshold(atm,h);
  std::cout << "Sim: h=" << h/100000 << ", e=" << e << ", eth=" << eth << '\n';
  int ilayer = det.getLayerNumber(h);
  unsigned lognprint=0;
  for(unsigned isim=3*nsim;isim;isim/=10)lognprint++;
  unsigned nprint=1;
  for(unsigned logiprint=2;logiprint<lognprint;logiprint++)nprint*=10;
  for(unsigned isim=0;isim<nsim;isim++)
    {
      egs5->shower(e, h/uz*ux, h/uz*uy, h, ux, uy, uz, ilayer, q);
      if(isim%nprint == (nprint-1))
	std::cout << "- nshower: " << isim+1 << '\n';
    }

  s->writeScalar("e",e);
  s->writeScalar("eth",eth);
  s->writeScalar("h",h);
  det.writeResults(s);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  // --------------------------------------------------------------------------
  // PROCESS OPTIONS
  // --------------------------------------------------------------------------

  bool print_usage = false;
  options.findBoolValue("h", print_usage, true, "Print this message.");
  bool print_usage2 = false;
  options.findBoolValue("help",  print_usage2, true, "Print this message.");
  print_usage |= print_usage2;
  std::string results_file("results.h5"); 
  options.findWithValue("o", results_file, "Output file name."); 

  options.addCatagory("primary", "Primary particle");
  double theta = 0; // degrees
  options.findWithValue("zn", theta, 
			"Set zenith angle of initial particle [deg].",
			"primary");

  double phi = 0; // degrees
  options.findWithValue("phi", phi, "Set azimuth angle (East to North) of "
			"initial particle [deg].", "primary");
  int q = -1; 
  options.findWithValue("q", q, "Set the initial particle type: -1=electron, "
			"0=photon or +1=positron.", "primary");
  double e = 40; // MeV
  options.findWithValue("e", e, "Set energy of initial particle [MeV].",
			"primary");
  
  double h = 10; // kilometers
  options.findWithValue("z", h, "Set height of initial particle [km].", 
			"primary");

  options.addCatagory("mode", "Running mode");
  std::string mode = "normal";

  triple<double,double,double> escan(20,2.5,50);
  if(options.findWithValue("escan", escan, "Enable energy scanning mode and "
			   "set parameters of scan. Values are Elo/dE/Ehi "
			   "[MeV]. If the parameters are negative then they "
			   "are taken to be relative to the Cherenkov "
			   "threshold.", "mode") != VSOptions::FS_NOT_FOUND)
    mode = "escan";

  triple<double,double,double> hscan(20,2.5,50);
  if(options.findWithValue("zscan", hscan, "Enable eltitude scanning mode and "
			   "set parameters of scan. Values are zlo/dz/zhi "
			   "[km].", "mode") != VSOptions::FS_NOT_FOUND)
    mode = "zscan";

  options.addCatagory("env", "Environment");
  std::string atm_file("atmprof6.dat");
  options.findWithValue("atm", atm_file, "File name of tabulated atmospheric "
			"model.", "env");
  bool enable_bfield = false;
  triple<double,double,double> bvec(-2831.9,11719.0,-25688.5);
  options.findBoolValue("bfield", enable_bfield, true, 
			"Enable magnetic field.", "env");
  options.findWithValue("bvec", bvec, 
			"Set the magnetic field vector EW/NS/UD [nT].", "env");

  unsigned nsim = 1000;
  options.findWithValue("n", nsim, "Number of simulations.");
 
  // --------------------------------------------------------------------------
  // FINISH OPTIONS PROCESSING
  // --------------------------------------------------------------------------

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

    if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  if(argc > 0)
    {
      std::cerr << progname << ": unrecognised parameters." << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }
  
  // --------------------------------------------------------------------------
  // GET ON WITH SIMULATION
  // --------------------------------------------------------------------------

  h *= 100000; // to CM
  hscan.first  *= 100000; // to CM
  hscan.second *= 100000; // to CM
  hscan.third  *= 100000; // to CM
  q = std::max(std::min(q,1),-1);
  bvec.first  *= 1e-5;
  bvec.second *= 1e-5;
  bvec.third  *= 1e-5;

  try
    {
      LayeredAtmosphere atm(atm_file);

      if(mode == "escan")
	{
	  double eth = cherenkovThreshold(atm,h);
	  if(escan.first < 0)escan.first *= -eth;
	  if(escan.second < 0)escan.second *= -eth;
	  if(escan.third < 0)escan.third *= -eth;
	}

      SimpleRNG* base_rng = EGS5RanluxSimpleRNG::instance();;

      //RandomNumbers rng_base(RandomNumbers::defaultFilename());
      //SimpleRNGAdapter rng(&rng_base);
      //base_rng = &rng;

      ConstantBField const_bfield(bvec.first, bvec.second, bvec.third);
      BField* bfield = enable_bfield ? &const_bfield : 0;
      double zbottom = 0;
      double ztop = HUGE_VAL;
      double zres = 10000;
      double yieldres = 1.0;
      ACDYield det(atm, 100, std::max(e,1000.0), bfield, yieldres, 
		   zres, zbottom, ztop);
      EGS5System* egs5 = EGS5System::instance(&det, base_rng);
      egs5->initializeEGS5();

      VSOctaveH5Writer writer(results_file, true);
      VSOctaveH5WriterStruct* settings = writer.writeStruct("settings");
      settings->writeScalar("theta",theta);
      settings->writeScalar("phi",phi);
      settings->writeScalar("q",q);
      settings->writeString("atm_file",atm_file);
      settings->writeScalar("bfield",enable_bfield);
      settings->writeScalar("bx",bvec.first);
      settings->writeScalar("by",bvec.second);
      settings->writeScalar("bz",bvec.third);
      settings->writeString("mode",mode);
      if(mode == "escan")
	{
	  settings->writeScalar("escan_Elo",escan.first);
	  settings->writeScalar("escan_Ehi",escan.third);
	  settings->writeScalar("escan_dE",escan.second);
	}
      else if(mode == "zscan")
	{
	  settings->writeScalar("zscan_zlo",hscan.first);
	  settings->writeScalar("zscan_zhi",hscan.third);
	  settings->writeScalar("zscan_dz",hscan.second);
	}
      delete settings;

      unsigned nloop = 0;
      double de = 0;
      double dh = 0;
      if(mode == "escan")
	{
	  e = escan.first;
	  de = escan.second;
	  if(de>0)
	    while(e<=escan.third)nloop++,e+=de;
	  else if(de<0)
	    while(e>=escan.third)nloop++,e+=de;
	  e = escan.first;
	}
      else if(mode == "zscan")
	{
	  h = hscan.first;
	  dh = hscan.second;
	  if(dh>0)
	    while(h<=hscan.third)nloop++,h+=dh;
	  else if(de<0)
	    while(h>=hscan.third)nloop++,h+=dh;
	  h = hscan.first;
	}

      double ux = -std::sin(theta/180*M_PI)*std::cos(phi/180.0*M_PI);
      double uy = -std::sin(theta/180*M_PI)*std::sin(phi/180.0*M_PI);
      double uz = -std::cos(theta/180*M_PI);
      
      if(nloop)
	{
	  std::vector<double> ve(nloop);
	  std::vector<double> veth(nloop);
	  std::vector<double> verel(nloop);
	  std::vector<double> vh(nloop);
	  std::vector<double> vymean(nloop);
	  std::vector<double> vyrms(nloop);

	  VSOctaveH5WriterCellVector* cv = 
	    writer.writeCellVector("scan_results", nloop);
	  for(unsigned iloop = 0; iloop<nloop; iloop++)
	    {
	      det.reset();
	      VSOctaveH5WriterStruct* s = cv->writeStruct(iloop);
	      runSim(h, e, ux, uy, uz, q, atm, egs5, det, nsim, s);
	      double eth = cherenkovThreshold(atm, h);
	      ve[iloop] = e;
	      veth[iloop] = eth;
	      verel[iloop] = e/eth;
	      vh[iloop] = h;
	      det.yieldStat(vymean[iloop], vyrms[iloop]);
	      delete s;
	      h += dh;
	      e += de;
	    }	  
	  delete cv;
	  VSOctaveH5WriterStruct* s = writer.writeStruct("scan_summary");
	  s->writeVector("e",ve);
	  s->writeVector("eth",veth);
	  s->writeVector("erel",verel);
	  s->writeVector("h",vh);
	  s->writeVector("yield_mean",vymean);
	  s->writeVector("yield_std",vyrms);
	  delete s;
	}
      else
	{
	  runSim(h, e, ux, uy, uz, q, atm, egs5, det, nsim, &writer);
	}
    }
  catch(std::string& s)
    {
      std::cerr << s << '\n';
    }
  catch(const char* s)
    {
      std::cerr << s << '\n';
    }
  std::cerr << "Fin..." << std::endl;
}

