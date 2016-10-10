// egs5_system.h - C header interface to EGS5. Structures are used to 
// represent (some of the) EGS5 common blocks.
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: egs5_system.h 4530 2012-09-27 08:42:49Z sfegan $

#ifndef EGS5_SYSTEM_H
#define EGS5_SYSTEM_H

#ifdef  __cplusplus
extern "C" {
#endif

// egrep '^[^!].*parameter' include/egs5_h.f | sed -e 's/^.*parameter (/#define /;s/=/ /; s/).*//' > egs5_params.h
#include "egs5_params.h"

// ! Cutoff energies & vacuum transport distance
struct EGS5_bounds
{
  double ecut[MXREG];     //          ! Charged particle cutoff energy (total)
  double pcut[MXREG];     //                    ! Photon cutoff energy (total)
  double vacdst;          //                       ! Vacuum transport distance
};

// ! Electron-photon control variables
struct EGS5_epcont    
{
  double edep;          //                                  ! Energy deposited
  double tstep;         //                      ! Distance to next interaction
  double ustep;         // ! User (straight line) step requested (and granted)
  double tmstep;        //          ! Total multiple scattering hinge distance
  double rhof;          //              ! Value density correction (default=1)
  double eold;          //  ! Charged particle (total) energy at start of step
  double enew;          //    ! Charged particle (total) energy at end of step
  double eke;           //                ! Kinetic energy of charged particle
  double elke;          //                          ! Natural logarithm of EKE
  double beta2;         //    ! Beta (v/c) squared of present charged particle
  double gle;           //                ! Natural logarithm of photon energy
  int idisc;            // ! User-discard flag: 0=no, >0=immediately, <0=later
  int irold;            //                          ! Index of previous region
  int irnew;            //                              ! Index of new regions
  int iausfl[MXAUS];    //              ! Flags for turning on calls to AUSGAB
};

// ! Names and data for media currently being used
struct EGS5_media
{
  double rlcm[MXMED];
  double rldu[MXMED];
  double rhom[MXMED];
  double charD[MXMED];
  int iraylm[MXMED];
  int incohm[MXMED];
  int iprofm[MXMED];
  int impacm[MXMED];
  int useGSD[MXMED];
  int media[MXMED][24];
  int nmed;
};

// ! Miscellaneous COMMON
struct EGS5_misc
{
  double rhor[MXREG];
  double k1Hscl[MXREG];
  double k1Lscl[MXREG];
  double ectrng[MXREG];
  double pctrng[MXREG];
  double dunit;
  int med[MXREG];
  int iraylr[MXREG];
  int lpolar[MXREG];
  int incohr[MXREG];
  int iprofr[MXREG];
  int impacr[MXREG];
  int nomsct[MXREG];      // ! For turning off mult. scattering in each region
  int kmpi;
  int kmpo;
  int nreg;               //               ! Number of regions in this problem
};

// ! Information kept about current particles
struct EGS5_stack
{
  double e[MXSTACK];      //  ! Total energy of particle (including rest mass)
  double x[MXSTACK];      //                                      ! X-position
  double y[MXSTACK];      //                                      ! Y-position
  double z[MXSTACK];      //                                      ! Z-position
  double u[MXSTACK];      //                         ! X-axis direction cosine
  double v[MXSTACK];      //                         ! Y-axis direction cosine
  double w[MXSTACK];      //                         ! Z-axis direction cosine
  double uf[MXSTACK];     //      ! Electric field vectors of polarized photon
  double vf[MXSTACK];
  double wf[MXSTACK];
  double dnear[MXSTACK];  //          ! Estimated distance to nearest boundary
  double wt[MXSTACK];     //                                 ! Particle weight
  double k1step[MXSTACK]; //                        ! Scat stren to next hinge
  double k1rsd[MXSTACK];  //            ! Scat stren from hinge to end of step
  double k1init[MXSTACK]; //                  ! Scat of prev hinge end of step
  double time[MXSTACK];
  double deinitial;
  double deresid;
  double denstep;
  int iq[MXSTACK];        //     ! Particle charge, -1[e-]; 0[photons]; +1(e+)
  int ir[MXSTACK];        //                                   ! Region number
  int latch[MXSTACK];     //                               ! Latching variable
  int np;                 //                             ! Stack pointer index
  int latchi;             //               ! Initialization for latch variable
};

// ! Some RANLUX settings
struct EGS5_rluxdat
{
  float twom12;
  float twom24;
  int ndskip;
  int luxlev;
  int nskip[5];
  int inseed;
  int kount;
  int mkount;
  int isdext[25];
  int rluxset;
};

// ! Some heavily used variables
struct EGS5_useful
{
  double rm;            //                                ! Electron rest mass
  int medium;           //            ! Index of current medium (0 for vacuum)
  int medold;           //                          ! Index of previous medium
  int iblobe;           //  ! Flag, photon below EBINDA after PE (1=yes, 0=no)
};

// ! User-Step-Controls
struct EGS5_usersc
{
  double estepr[MXREG];   // ! estepe multiplier for each region (if non zero)
  double esave[MXREG];    //         ! Upper limit on electron range rejection
  double emaxe;           //               ! maximum kinetic energy in problem
};

// ! Photon transport data
struct EGS5_photin
{
  double ebinsa[100];             // Average K-edge binding energy
  double ge0[MXMED];
  double ge1[MXMED];
  double gmfp0[MXMED][MXGE];      // Gamma Mean Free Path fit
  double gmfp1[MXMED][MXGE];      //  coefficients
  double gbr10[MXMED][MXGE];      // Gamma Branching Ratio fit coefficients
  double gbr11[MXMED][MXGE];      //  (pair/(pair+compton+photo)
  double gbr20[MXMED][MXGE];      // Gamma Branching Ratio fit coefficients
  double gbr21[MXMED][MXGE];      //  (pair+compton)/(pair+compton+photo)
  double rco0[MXMED];
  double rco1[MXMED];
  double rsct0[MXMED][MXRAYFF];
  double rsct1[MXMED][MXRAYFF];
  double cohe0[MXMED][MXGE];
  double cohe1[MXMED][MXGE];
  int mpgem[MXMED][MXSGE];
  int ngr[MXMED];
};

// ! Electron transport input
struct EGS5_elecin
{
  double ekelim;
  double eke0[MXMED];
  double eke1[MXMED];
  double cmfp0[MXMED],cmfp1[MXMED];
  double xr0[MXMED],teff0[MXMED],blcc[MXMED],xcc[MXMED];
  double picmp0[MXMED][MXCMFP],picmp1[MXMED][MXCMFP];
  double eicmp0[MXMED][MXCMFP],eicmp1[MXMED][MXCMFP];
  double esig0[MXMED][MXEKE],esig1[MXMED][MXEKE];
  double psig0[MXMED][MXEKE],psig1[MXMED][MXEKE];
  double ededx0[MXMED][MXEKE],ededx1[MXMED][MXEKE];
  double pdedx0[MXMED][MXEKE],pdedx1[MXMED][MXEKE];
  double erang0[MXMED][MXEKE],erang1[MXMED][MXEKE];
  double prang0[MXMED][MXEKE],prang1[MXMED][MXEKE];
  double estep0[MXMED][MXEKE],estep1[MXMED][MXEKE];
  double ebr10[MXMED][MXEKE],ebr11[MXMED][MXEKE];
  double pbr10[MXMED][MXEKE],pbr11[MXMED][MXEKE];
  double pbr20[MXMED][MXEKE],pbr21[MXMED][MXEKE];
  double tmxs0[MXMED][MXEKE],tmxs1[MXMED][MXEKE];
  double cmfpe0[MXMED][MXLEKE],cmfpe1[MXMED][MXLEKE];
  double cmfpp0[MXMED][MXLEKE],cmfpp1[MXMED][MXLEKE];
  double cxc2e0[MXMED][MXLEKE],cxc2e1[MXMED][MXLEKE];
  double cxc2p0[MXMED][MXLEKE],cxc2p1[MXMED][MXLEKE];
  double clxae0[MXMED][MXLEKE],clxae1[MXMED][MXLEKE];
  double clxap0[MXMED][MXLEKE],clxap1[MXMED][MXLEKE];
  double thr0[MXRNTH][MXBLC],thr1[MXRNTH][MXBLC],thr2[MXRNTH][MXBLC];
  double thri0[MXRNTHI][MXBLC],thri1[MXRNTHI][MXBLC],thri2[MXRNTHI][MXBLC];
  double fstep[MSSTEPS],fsqr[MSSTEPS];
  double vert1[MXVRT1],vert2[MSSTEPS][MXVRT2];
  double blc0,blc1,rthr0,rthr1,rthri0,rthri1;
  int icomp;
  int mpeem[MXMED][MXSEKE],msmap[MXJREFF];
  int msteps,jrmax,mxv1,mxv2,nblc,nrnth,nrnthi;
  int iunrst[MXMED],epstfl[MXMED],iaprim[MXMED];
};

// ! Bremsstrahlung and pair production data
struct EGS5_brempr
{
  unsigned asym[2][MXEL][MXMED];
  double dl1[MXMED][6];
  double dl2[MXMED][6];
  double dl3[MXMED][6];
  double dl4[MXMED][6];
  double dl5[MXMED][6];
  double dl6[MXMED][6];
  double alphi[MXMED][2];
  double bpar[MXMED][2];
  double delpos[MXMED][2];
  double wa[MXEL][MXMED];
  double pz[MXEL][MXMED];
  double zelem[MXEL][MXMED];
  double rhoz[MXEL][MXMED];
  double pwr2i[MXPWR2I];
  double delcm[MXMED];
  double zbrang[MXMED];
  double fbrspl;
  int nne[MXMED];
  int ibrdst;
  int iprdst;
  int ibrspl;
  int nbrspl;
};

extern void block_set_(void);
extern void counters_out_(int*);
extern void hatch_(void);
extern void pegs5_(void);
extern void rluxinit_(void);
extern void shower_(int* iqi, double* ei,
		    double* xi,double* yi,double* zi,
		    double* ui,double* vi,double* wi, 
		    int *iri, double *wti);
extern void hardx_(int*,double*,int*,double*,double*);
extern void randomset_(double*,int*);
extern void origrandomset_(double*);
extern void ncallrandomset_(int*);

// User defined functions that must be specified
extern void ausgab_(int* iarg);
extern void howfar_(void);

#ifdef  __cplusplus
} // extern "C"
#endif

#endif // defined EGS5_SYSTEM_H
