// image.cpp - General test program for development
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: test.cpp 4802 2012-11-22 14:58:23Z sfegan $

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

#include <VSOptions.hpp>
#include <VSOctaveH5Writer.hpp>

#include "EGS5LayeredDetector.hpp"
#include "EGS5AtmosphericDetector.hpp"
#include "DetectorEfficiency.hpp"
#include "Atmosphere.hpp"

#include "VSSimpleHist.hpp"
#include "RandomNumbers.hpp"
#include "MultiRNG.hpp"

using namespace VERITAS;

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
  std::string image_file("doodles/image.h5"); 
  options.findWithValue("o", image_file, "Output file name."); 

  // Primary ------------------------------------------------------------------

  options.addCatagory("primary", "Primary particle");
  double theta = 0; // degrees
  options.findWithValue("zn", theta, 
			"Set zenith angle of initial particle [deg].",
			"primary");

  double az = 0; // degrees
  options.findWithValue("az", az, "Set azimuth angle (North to East) of "
			"initial particle [deg].", "primary");
  int q = 0; 
  options.findWithValue("q", q, "Set the initial particle type: -1=electron, "
			"0=photon or +1=positron.", "primary");
  double e = 100; // GeV
  options.findWithValue("e", e, "Set energy of initial particle [GeV].",
			"primary");
  
  double h = HUGE_VAL; // kilometers
  options.findWithValue("z", h, "Set height of initial particle [km].", 
			"primary");

  // Mode ---------------------------------------------------------------------

  //options.addCatagory("mode", "Running mode");

  // Environment --------------------------------------------------------------

  options.addCatagory("env", "Environment");
  std::string atm_file("Parameters/atmprof6.dat");
  options.findWithValue("atm", atm_file, "File name of tabulated atmospheric "
			"model.", "env");
  unsigned nlayer = 100;
  options.findWithValue("nlayer", nlayer, "Number of layeres to use in"
			"layered atmosphere.", "env");
  unsigned nmedia = 1;
  options.findWithValue("nmedia", nmedia, "Number of air media to use in "
			"layered atmosphere.", "env");
  bool enable_bfield = false;
  triple<double,double,double> bvec(-2831.9,11719.0,-25688.5);
  options.findBoolValue("bfield", enable_bfield, true, 
			"Enable magnetic field.", "env");
  options.findWithValue("bvec", bvec, 
			"Set the magnetic field vector EW/NS/UD [nT].", "env");

  // Instrument ---------------------------------------------------------------

  options.addCatagory("instrument", "Instrument");
  std::string mirr_eff_file("Parameters/corsika_mirreff.dat");
  std::string quant_eff_file("Parameters/corsika_quanteff.dat");
  std::string atm_abs_file("Parameters/corsika_atmabs.dat");  

  options.findWithValue("mirr_eff", mirr_eff_file, "Wavelength dependent "
			"mirror reflectivity file.", "instrument");
  options.findWithValue("quant_eff", quant_eff_file, "Wavelength dependent "
			"quantum efficienty file.", "instrument");
  options.findWithValue("atm_abs", atm_abs_file, "Wavelength and altitude "
			"dependent atmospheric absorption file.",
			"instrument");

  double z0 = 1.5;
  options.findWithValue("z0", z0, "Altitude of the telescope array [km].",
			"instrument");

  typedef triple<double,double,double> ss_type;
  std::vector<ss_type> scope_spec;
  scope_spec.push_back(ss_type( 60, 60, 12));
  scope_spec.push_back(ss_type(-60, 60, 12));
  scope_spec.push_back(ss_type( 60,-60, 12));
  scope_spec.push_back(ss_type(-60,-60, 12));
  scope_spec.push_back(ss_type(  0,  0, 28));
  options.findWithValue("scopes", scope_spec, "Positions and diameters of "
			"telescopes in array. Vector of triples of X,Y,D, "
			"where X is E-W, Y is N-S [m].", "instrument");

  double res = 0.02;
  options.findWithValue("res", res, "Resolution of the telescope [deg]",
			"instrument");
  double fov = 5.0;
  options.findWithValue("fov", fov, "Field of view of the telescope [deg]",
			"instrument");

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

  if(h != HUGE_VAL)h *= 100000; // from km to cm
  z0 *= 100000; // from km to cm
  q = std::max(std::min(q,1),-1);
  bvec.first  *= 1e-5; // from T to G
  bvec.second *= 1e-5; // from T to G
  bvec.third  *= 1e-5; // from T to G
  theta *= M_PI/180.0;
  az *= M_PI/180.0;
  e *= 1000.0;

  try
    {
      TelescopeEfficiency eff;
      eff.scaleEffFromFile(mirr_eff_file);
      eff.scaleEffFromFile(quant_eff_file);
      AtmosphericAbsorption atmabs(atm_abs_file);

      SimpleRNG* base_rng = EGS5RanluxSimpleRNG::instance();

      //RandomNumbers rng_base(RandomNumbers::defaultFilename());
      //SimpleRNGAdapter rng(&rng_base);
      //base_rng = &rng;

      PredefinedDeviateSimpleRNG pd_rng(base_rng);
      pd_rng.setPD(0,0.3);

      ConstantDeviateSimpleRNG c_rng(0.5);

      LayeredAtmosphere atm(atm_file);
      double zbottom = z0;
      double ztop = std::min(h, atm.topOfAtmosphere());

      ConstantBField const_bfield(bvec.first, bvec.second, bvec.third);
      BField* bfield = enable_bfield ? &const_bfield : 0;

      EGS5SimpleIACTArray det(atm, nlayer, 1E6, bfield, zbottom, ztop, nmedia);
      std::vector<Layer> layers;
      det.getLayers(layers);

      for(std::vector<ss_type>::const_iterator iscope = scope_spec.begin();
	  iscope != scope_spec.end(); iscope++)
	{
	  EGS5SimpleIACTArray::ImagingScope s; 
	  s.x.set(iscope->first * 100, iscope->second * 100, z0);
	  s.r   = iscope->third * 100 * 0.5;
	  s.res = res;
	  s.fov = fov;
	  det.addScope(s);
	}

      double ctheta = std::cos(theta);
      double stheta = std::sin(theta);

      for(std::vector<EGS5SimpleIACTArray::ImagingScope>::iterator
	    iscope = det.scopes().begin(); iscope != det.scopes().end();
	  iscope++)
	{
	  iscope->setZnAz(theta, az);
	  iscope->yield = atmabs.integrateYield(iscope->x.z(), ctheta, eff);
	}

      EGS5System* egs5 = EGS5System::instance(&det, base_rng);
      //egs5->setRNG(35, &c_rng); // Number of MFP for hard step in e+/e-
      //egs5->setRNG(86, &pd_rng); // Number of MFP for photon      
      egs5->initializeEGS5();
      
      for(unsigned i=0;i<nsim;i++)
	{
	  pd_rng.resetSession();
	  egs5->shower(e, 6000, 0, layers.back().zt,
		       stheta, 0.0, -ctheta,
		       layers.size()+1, q);
	  std::cout << i+1 << '\n';
	}

      unsigned nscope = det.images().size();
      VSOctaveH5Writer writer(image_file,true);
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
