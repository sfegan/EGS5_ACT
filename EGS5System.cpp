// EGS5System.cpp - Wrapper around EGS5 FORTRAN subroutines
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5System.cpp 5698 2013-10-02 13:10:45Z sfegan $

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fcntl.h>
#include <time.h>

#include "EGS5System.hpp"

#define stringify(c) xstringify(c)
#define xstringify(s) #s

// ****************************************************************************
//
// EGS5UserInterface
//
// ****************************************************************************

EGS5UserInterface::~EGS5UserInterface()
{
  // nothing to see here
}

void EGS5UserInterface::egs5InitializeStarting()
{
  // nothing to see here
}

void EGS5UserInterface::egs5InitializeCompleted()
{
  // nothing to see here
}

bool EGS5UserInterface::EGS5UserInterface::shouldRunPEGS5()
{
  return true;
}

void EGS5UserInterface::writePEGS5InputFile(const std::string& filename)
{
  // nothing to see here
}

void EGS5UserInterface::showerStarting(double e,
				       double x, double y, double z,
				       double u, double v, double w,
				       int iregion, int q, double weight)
{
  // nothing to see here
}

void EGS5UserInterface::showerCompleted()
{
  // nothing to see here
}

// ****************************************************************************
//
// EGS5UIDelegator
//
// ****************************************************************************

EGS5UIDelegator::EGS5UIDelegator(): EGS5UserInterface(), m_delegatees()
{
  // nothing to see here
}

EGS5UIDelegator::~EGS5UIDelegator()
{
  // nothing to see here
}

unsigned EGS5UIDelegator::
appendDelegatee(EGS5UserInterface* ui, bool delegate_all)
{
  Delegatee d;
  d.ui = ui;
  d.delegateGetEMax = delegate_all;
  d.delegateGetMedia = delegate_all;
  d.delegateGetRegions = delegate_all;
  d.delegateAusgab = delegate_all;
  d.delegateHowfar = delegate_all;
  d.delegateEGS5InitializeStarting = delegate_all;
  d.delegateEGS5InitializeCompleted = delegate_all;
  d.delegateShouldRunPEGS5 = delegate_all;
  d.delegateWritePEGS5InputFile = delegate_all;
  d.delegateShowerStarting = delegate_all;
  d.delegateShowerCompleted = delegate_all;
  d.stopAfterDiscardInHowfar = true;
  m_delegatees.push_back(d);
  return m_delegatees.size()-1;
}

void EGS5UIDelegator::setDelegateAll(unsigned id, bool enable)
{
  setDelegateGetEMax(id, enable);
  setDelegateGetMedia(id, enable);
  setDelegateGetRegions(id, enable);
  setDelegateAusgab(id, enable);
  setDelegateHowfar(id, enable);
  setDelegateEGS5InitializeStarting(id, enable);
  setDelegateEGS5InitializeCompleted(id, enable);
  setDelegateShouldRunPEGS5(id, enable);
  setDelegateWritePEGS5InputFile(id, enable);
  setDelegateShowerStarting(id, enable);
  setDelegateShowerCompleted(id, enable);
}

void EGS5UIDelegator::setDelegateGetEMax(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateGetEMax = enable;
}

void EGS5UIDelegator::setDelegateGetMedia(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateGetMedia = enable;
}

void EGS5UIDelegator::setDelegateGetRegions(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateGetRegions = enable;
}

void EGS5UIDelegator::setDelegateAusgab(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateAusgab = enable;
}

void EGS5UIDelegator::setDelegateHowfar(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateHowfar = enable;
}

void EGS5UIDelegator::setDelegateEGS5InitializeStarting(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateEGS5InitializeStarting = enable;
}

void EGS5UIDelegator::setDelegateEGS5InitializeCompleted(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateEGS5InitializeCompleted = enable;
}

void EGS5UIDelegator::setDelegateShouldRunPEGS5(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateShouldRunPEGS5 = enable;
}

void EGS5UIDelegator::setDelegateWritePEGS5InputFile(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateWritePEGS5InputFile = enable;
}

void EGS5UIDelegator::setDelegateShowerStarting(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateShowerStarting = enable;
}

void EGS5UIDelegator::setDelegateShowerCompleted(unsigned id, bool enable)
{
  m_delegatees.at(id).delegateShowerCompleted = enable;
}

void EGS5UIDelegator::setStopAfterDiscardInHowfar(unsigned id, bool enable)
{
  m_delegatees.at(id).stopAfterDiscardInHowfar = enable;
}

double EGS5UIDelegator::getEMax()
{
  double emax = 0;
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateGetEMax)
      emax = std::max(emax, id->ui->getEMax());
  return emax;
}

void EGS5UIDelegator::getMedia(std::vector<Media>& media)
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateGetMedia)
      id->ui->getMedia(media);
}

void EGS5UIDelegator::getRegions(std::vector<Region>& regions)
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateGetRegions)
      id->ui->getRegions(regions);
}

void EGS5UIDelegator::ausgab(int& iarg)
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateAusgab)
      id->ui->ausgab(iarg);
}

void EGS5UIDelegator::howfar(void)
{
  extern EGS5_epcont epcont_;
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    {
      if(id->delegateHowfar)id->ui->howfar();
      if(epcont_.idisc!=0 && id->stopAfterDiscardInHowfar)break;
    }
}

void EGS5UIDelegator::egs5InitializeStarting()
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateEGS5InitializeStarting)
      id->ui->egs5InitializeStarting();
}

void EGS5UIDelegator::egs5InitializeCompleted()
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateEGS5InitializeCompleted)
      id->ui->egs5InitializeCompleted();
}

bool EGS5UIDelegator::shouldRunPEGS5()
{
  bool run = false;
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateShouldRunPEGS5)
      run |= id->ui->shouldRunPEGS5();
  return run;
}

void EGS5UIDelegator::writePEGS5InputFile(const std::string& filename)
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateWritePEGS5InputFile)
      id->ui->writePEGS5InputFile(filename);
}


void EGS5UIDelegator::showerStarting(double e,
				     double x, double y, double z,
				     double u, double v, double w,
				     int iregion, int q, double weight)
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateShowerStarting)
      id->ui->showerStarting(e,x,y,z,u,v,w,iregion,q,weight);
}

void EGS5UIDelegator::showerCompleted()
{
  for(std::vector<Delegatee>::iterator id=m_delegatees.begin(); 
      id != m_delegatees.end(); id++)
    if(id->delegateShowerCompleted)
      id->ui->showerCompleted();
}

// ****************************************************************************
//
// EGS5RanluxSimpleRNG
//
// ****************************************************************************

static std::auto_ptr<EGS5RanluxSimpleRNG> s_egs_ranlux_singleton;

EGS5RanluxSimpleRNG::~EGS5RanluxSimpleRNG()
{
  // nothing to see here
}

double EGS5RanluxSimpleRNG::uniform()
{
  double d = 0;
  origrandomset_(&d);
  return d;
}

void EGS5RanluxSimpleRNG::setSeed(int seed)
{
  // Initialize the RNG
  seed &= 0x7FFFFFFF;
  if(seed==0)
    {
      int fd = open("/dev/urandom",O_RDONLY);      
      if(fd>=0)
	{
	  while(seed==0)
	    {
	      read(fd,&seed,sizeof(seed));
	      seed &= 0x7FFFFFFF;
	    }
	  close(fd);
	}
      else
	{
	  time_t t;
	  time(&t);
	  seed = t&0x7FFFFFFF;
	}
    }
  std::cerr << "Seed: " << seed << std::endl;
  extern struct EGS5_rluxdat rluxdat_;
  rluxdat_.inseed=seed;
  rluxdat_.luxlev=1;
  rluxinit_();  
}

EGS5RanluxSimpleRNG* EGS5RanluxSimpleRNG::instance(int seed)
{
  EGS5RanluxSimpleRNG* s = s_egs_ranlux_singleton.get();
  if(s==0)
    {
      s = new EGS5RanluxSimpleRNG();
      s->setSeed(seed);
      s_egs_ranlux_singleton.reset(s);
    }
  return s;
}

// ****************************************************************************
//
// EGS5System
//
// ****************************************************************************

static std::auto_ptr<EGS5System> s_egs_system_singleton;

EGS5System* EGS5System::instance(EGS5UserInterface* _ui, SimpleRNG* _rng)
{
  EGS5System* s = s_egs_system_singleton.get();
  if(s==0)
    {
      s = new EGS5System();
      s_egs_system_singleton.reset(s);
    }
  if(_ui)s->setUI(_ui);
  if(_rng)s->setAllRNG(_rng);
  return s;
}

void EGS5System::setUI(EGS5UserInterface* _ui, bool force_reinitialize)
{
  m_ui = _ui;
  if(force_reinitialize)m_initialized = false;
}

void EGS5System::setRNG(unsigned icall, SimpleRNG* _rng)
{
  m_rng.at(icall) = _rng;
}

void EGS5System::setAllRNG(SimpleRNG* _rng)
{
  for(unsigned icall=0;icall<m_rng.size();icall++)m_rng[icall] = _rng;
}

void EGS5System::setEnableBremsPPThetaSampling(bool enable)
{
  if(m_initialized)
    std::cout << "Warning: enabling brems and PP theta sampling when already initialized" << std::endl;
  m_enable_bremspp_theta_sampling = enable;
}

void EGS5System::initializeEGS5()
{  
  if(!m_ui)throw std::string("EGS5System: no user interface specified");

  extern struct EGS5_usersc usersc_;

  // Let the UI know we are initializing EGS5
  m_ui->egs5InitializeStarting();

  block_set_();
   
  // Get the list of media from the UI and run PEGS5 if requested
  std::vector<EGS5UserInterface::Media> media;
  m_ui->getMedia(media);
  if(media.size() == 0)
    throw std::string("EGS5System: must supply at least one media");
  if(media.size() > MXMED)
    throw std::string("EGS5System: number of media exceeds compile-time "
		      "limit of " stringify(MXMED));

  // Get the list of media from the UI and run PEGS5 if requested
  extern struct EGS5_media media_;
  media_.nmed=media.size();
  for(int imed=0;imed<media_.nmed;imed++)
    {
      unsigned ichar=0;
      for(;ichar<24 && ichar<media[imed].name.size(); ichar++)
	media_.media[imed][ichar] = 0x20202000|int(media[imed].name[ichar]);
      for(;ichar<24; ichar++)
	media_.media[imed][ichar] = 0x20202020;
      media_.charD[imed] = media[imed].charD;
    }

  if(m_ui->shouldRunPEGS5())
    {
      m_ui->writePEGS5InputFile("pgs5job.pegs5inp");
      pegs5_();
    }

  // Get the list of regions from the UI
  std::vector<EGS5UserInterface::Region> regions;
  m_ui->getRegions(regions);
  if(regions.size() == 0)
    throw std::string("EGS5System: must supply at least one region");
  if(regions.size() > MXREG)
    throw std::string("EGS5System: number of regions exceeds compile-time "
		      "limit of " stringify(MXREG));

  extern struct EGS5_misc misc_;
  extern struct EGS5_bounds bounds_;
  misc_.nreg = regions.size();
  for(unsigned iregion=0; iregion<regions.size(); iregion++)
    {
      if((regions[iregion].imedia<=-2)
	 ||(regions[iregion].imedia>=int(media.size())))
	{
	  std::ostringstream s;
	  s << "EGS5System: illegal media number: " << regions[iregion].imedia
	    << " in region: " << iregion;
	  throw s.str();
	}

      misc_.med[iregion]    = regions[iregion].imedia+1;
      misc_.rhor[iregion]   = regions[iregion].rho;
      bounds_.ecut[iregion] = regions[iregion].ecut;
      bounds_.pcut[iregion] = regions[iregion].pcut;
    }

  // Set the RNG if necessary - use the EGS5 generator as a fallback 
  SimpleRNG* rng = 0;
  for(unsigned icall=0;icall<m_rng.size();icall++)
    if(m_rng[icall]==0)
      {
	if(rng==0)rng = EGS5RanluxSimpleRNG::instance();
	m_rng[icall] = rng;
      }
  
  // Get the maximum energy from the UI
  usersc_.emaxe = m_ui->getEMax();
  
  // Enable sampling of Bremstrahlung and Pair Production angular distributions
  if(m_enable_bremspp_theta_sampling)
    {
      extern EGS5_brempr brempr_;
      brempr_.ibrdst = 1;
      brempr_.iprdst = 2;
    }

  // Open the PEGS5 data file and "hatch" the EGS5 system
  std::string fn_pegs5dat("pgs5job.pegs5dat");
  std::string fn_dummy("egs5job.dummy");
  fort_open_old_(&misc_.kmpi,fn_pegs5dat.c_str(),fn_pegs5dat.size());
  fort_open_unknown_(&misc_.kmpo,fn_dummy.c_str(),fn_dummy.size());
  hatch_();
  fort_close_(&misc_.kmpi);
  fort_close_(&misc_.kmpo);

  m_initialized = true;

  // Let the UI know we are finished inigtializing EGS5
  m_ui->egs5InitializeCompleted();
}

void EGS5System::shower(double e,
			double x, double y, double z,
			double u, double v, double w,
			int iregion, int q, double weight)
{			
  if(!m_initialized)
    throw std::string("EGS5System: must be initialized before calling shower");
  m_ui->showerStarting(e,x,y,z,u,v,w,iregion,q,weight);
  iregion++;
  shower_(&q,&e,&x,&y,&z,&u,&v,&w,&iregion,&weight);
  m_ui->showerCompleted();
}

EGS5System::EGS5System(): 
  m_rng(), m_ui(), m_enable_bremspp_theta_sampling(true), m_initialized(false)
{
  int i0 = 0;
  counters_out_(&i0);
  int i66 = 66;
  std::string fn("egs5job.log");
  fort_open_unknown_(&i66,fn.c_str(),fn.size());  
  int ncall;
  ncallrandomset_(&ncall);
  m_rng.resize(ncall);
}

EGS5System::~EGS5System()
{
  int i1 = 1;
  int i66 = 66;
  int i99 = 99;
  std::string fn("egs5job.counters");
  fort_open_unknown_(&i99,fn.c_str(),fn.size());
  counters_out_(&i1);
  fort_close_(&i99);
  fort_close_(&i66);
}

// ****************************************************************************
//
// FORTRAN interface functions
//
// ****************************************************************************

void ausgab_(int* iarg)
{
  s_egs_system_singleton->ui()->ausgab(*iarg);
}

void howfar_(void)
{
  s_egs_system_singleton->ui()->howfar();
}

void randomset_(double* d, int* icall)
{
  *d = s_egs_system_singleton->rng(*icall)->uniform();
}
