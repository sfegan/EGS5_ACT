// EGS5System.hpp - Wrapper around EGS5 FORTRAN subroutines
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5System.hpp 5698 2013-10-02 13:10:45Z sfegan $

#ifndef EGS5SYSTEM_HPP
#define EGS5SYSTEM_HPP

#include <vector>
#include <string>
#include <memory>

#include <SimpleRNG.hpp>

#include "egs5_system.h"
#include "egs5_futils.h"

#define E5S_RCURV_CM_PER_MEV_OVER_GAUSS 3335.64095198152

// ****************************************************************************
// ****************************************************************************
// ***                                                                      ***
// *** EGS5 USER INTERFACE                                                  ***
// ***                                                                      ***
// ****************************************************************************
// ****************************************************************************

struct EGS5Media
{
  EGS5Media(const std::string& _name = std::string(), double _charD = 0): 
    name(_name), charD(_charD) { }
  std::string name;
  double charD;
};
  
struct EGS5Region
{
  EGS5Region(int _imedia = 0, double _rho = 0, double _ecut = 0, 
	     double _pcut = 0): 
    imedia(_imedia), rho(_rho), ecut(_ecut), pcut(_pcut) { }
  int imedia;   // Index into media array - (-1 is vacuum)
  double rho;   // Override default media density
  double ecut;  // Override default electron cutoff
  double pcut;  // Override default photon cutoff
};

class EGS5UserInterface
{
public:
  typedef EGS5Media   Media;
  typedef EGS5Region  Region;

  virtual ~EGS5UserInterface();

  // **************************************************************************
  // REQUIRED FUNCTIONS - must be supplied by derived class
  // **************************************************************************

  // Return maximum energy
  virtual double getEMax() = 0;

  // Return list of media names
  virtual void getMedia(std::vector<Media>& media) = 0;

  // Return list of regions in problem
  virtual void getRegions(std::vector<Region>& regions) = 0;

  // Primary EGS5 interface
  virtual void ausgab(int& iarg) = 0;
  virtual void howfar(void) = 0;

  // **************************************************************************
  // OPTIONAL FUNCTIONS
  // **************************************************************************

  // EGS5System has been asked to initialize system - default do nothing
  virtual void egs5InitializeStarting();
  virtual void egs5InitializeCompleted();

  // Should EGS5System run PEGS5 - defaults to true
  virtual bool shouldRunPEGS5();

  // Write the PEGS5 input file - defaults to doing nothing
  virtual void writePEGS5InputFile(const std::string& filename);

  // Shower starting
  virtual void showerStarting(double e,
			      double x, double y, double z,
			      double u, double v, double w,
			      int iregion, int q, double weight);
  virtual void showerCompleted();
};

// ****************************************************************************
// ****************************************************************************
// ***                                                                      ***
// *** EGS5 USER INTERFACE DELEGATOR                                        ***
// ***                                                                      ***
// ****************************************************************************
// ****************************************************************************

class EGS5UIDelegator: public EGS5UserInterface
{
public:
  EGS5UIDelegator();
  virtual ~EGS5UIDelegator();

  unsigned appendDelegatee(EGS5UserInterface* ui, bool delegate_all=true);
  
  void setDelegateAll(unsigned id, bool enable=true);
  void setDelegateGetEMax(unsigned id, bool enable=true);
  void setDelegateGetMedia(unsigned id, bool enable=true);
  void setDelegateGetRegions(unsigned id, bool enable=true);
  void setDelegateAusgab(unsigned id, bool enable=true);
  void setDelegateHowfar(unsigned id, bool enable=true);
  void setDelegateEGS5InitializeStarting(unsigned id, bool enable=true);
  void setDelegateEGS5InitializeCompleted(unsigned id, bool enable=true);
  void setDelegateShouldRunPEGS5(unsigned id, bool enable=true);
  void setDelegateWritePEGS5InputFile(unsigned id, bool enable=true);
  void setDelegateShowerStarting(unsigned id, bool enable=true);
  void setDelegateShowerCompleted(unsigned id, bool enable=true);
  void setStopAfterDiscardInHowfar(unsigned id, bool enable=true);

  // EGS5 UI interface
  virtual double getEMax();
  virtual void getMedia(std::vector<Media>& media);
  virtual void getRegions(std::vector<Region>& regions);

  virtual void ausgab(int& iarg);
  virtual void howfar(void);

  virtual void egs5InitializeStarting();
  virtual void egs5InitializeCompleted();

  virtual bool shouldRunPEGS5();

  virtual void writePEGS5InputFile(const std::string& filename);

  virtual void showerStarting(double e,
			      double x, double y, double z,
			      double u, double v, double w,
			      int iregion, int q, double weight);
  virtual void showerCompleted();

private:
  class Delegatee
  {
  public:
    EGS5UserInterface* ui;
    bool delegateGetEMax;
    bool delegateGetMedia;
    bool delegateGetRegions;
    bool delegateAusgab;
    bool delegateHowfar;
    bool delegateEGS5InitializeStarting;
    bool delegateEGS5InitializeCompleted;
    bool delegateShouldRunPEGS5;
    bool delegateWritePEGS5InputFile;
    bool delegateShowerStarting;
    bool delegateShowerCompleted;
    bool stopAfterDiscardInHowfar;
  };

  std::vector<Delegatee> m_delegatees;
};

// ****************************************************************************
// ****************************************************************************
// ***                                                                      ***
// *** EGS5 RANLUX SIMPLE RNG ADAPTER                                       ***
// ***                                                                      ***
// ****************************************************************************
// ****************************************************************************

class EGS5RanluxSimpleRNG: public SimpleRNG
{
public:
  virtual ~EGS5RanluxSimpleRNG();
  virtual double uniform();
  void setSeed(int seed = 0);
  static EGS5RanluxSimpleRNG* instance(int seed = 0);
private:
  EGS5RanluxSimpleRNG(): SimpleRNG() { }
};

// ****************************************************************************
// ****************************************************************************
// ***                                                                      ***
// *** EGS5 SYSTEM                                                          ***
// ***                                                                      ***
// ****************************************************************************
// ****************************************************************************

class EGS5System
{
public:
  ~EGS5System();

  // Return singleton instance of the EGS5 system
  static EGS5System* instance(EGS5UserInterface* _ui=0, SimpleRNG* _rng=0);

  // Get and set the interface used by the EGS5 interface
  EGS5UserInterface* ui() const { return m_ui; }
  void setUI(EGS5UserInterface* _ui, bool force_reinitialize = true);  

  // Get and set the RNG
  unsigned nRNG() const { return m_rng.size(); }
  SimpleRNG* rng(unsigned icall) const { return m_rng[icall]; }
  void setAllRNG(SimpleRNG* _rng);
  void setRNG(unsigned icall, SimpleRNG* _rng);

  // Enable or disable "full" Bremstrahlung and PP theta distribution
  // Note: default is "enable"
  void setEnableBremsPPThetaSampling(bool enable = true);

  // Initialize the EGS5 system - define the regions and media etc
  void initializeEGS5();

  // Generate shower
  void shower(double e,
	      double x=0.0, double y=0.0, double z=0.0,
	      double u=0.0, double v=0.0, double w=1.0,
	      int iregion=0, int q=0, double weight=1.0);

protected:
  EGS5System();
  std::vector<SimpleRNG*>          m_rng;
  EGS5UserInterface*               m_ui;
  bool                             m_enable_bremspp_theta_sampling;
  bool                             m_initialized;
};

#endif // not defined EGS5SYSTEM_H
