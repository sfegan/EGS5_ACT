// Reference EGS5 test application - 1000 showers @ 100 GeV = 12 sec
// Stephen Fegan - 2012-10-25
// g++-mp-4.6 -O3 -o simplest_setup simplest_setup.cpp egs5_system.o egs5_futils.o -lgfortran

#include <cstdlib>
#include <cstring>

#include "egs5_system.h"
#include "egs5_futils.h"

int main()
{
  int izero = 0;
  int ione = 1;
  const char* media[] = { "AIR AT NTP" };

  extern struct EGS5_useful useful_;
  extern struct EGS5_media media_;
  extern struct EGS5_misc misc_;
  extern struct EGS5_rluxdat rluxdat_;
  extern struct EGS5_usersc usersc_;

  counters_out_(&izero);
  block_set_();

  media_.nmed=1;
  for(int imed=0;imed<media_.nmed;imed++)
    {
      unsigned ichar=0;
      while(ichar<24 && media[imed][ichar])
	media_.media[imed][ichar] = 0x20202000|int(media[imed][ichar]), ichar++;
      while(ichar<24)
	media_.media[imed][ichar] = 0x20202000|int(' '), ichar++;
    }
  media_.charD[0]=10000;
  pegs5_();

  misc_.nreg = 1;
  misc_.med[0] = 1;

  rluxdat_.inseed=1234567;
  rluxdat_.luxlev=1;
  rluxinit_();

  int iqi=0;
  double xi=0.0;
  double yi=0.0;
  double zi=0.0;
  double ui=0.0;
  double vi=0.0;
  double wi=1.0;
  int  iri=1;
  double wti=1.0;
  double ei=100000;

  if(iqi != 0)usersc_.emaxe = ei;
  else usersc_.emaxe = ei + useful_.rm;
  
  char fn_pegs5dat[] = "pgs5job.pegs5dat";
  char fn_dummy[] = "egs5job.dummy";
  fort_open_old_(&misc_.kmpi,fn_pegs5dat,strlen(fn_pegs5dat));
  fort_open_unknown_(&misc_.kmpo,fn_dummy,strlen(fn_dummy));
  hatch_();
  fort_close_(&misc_.kmpi);
  fort_close_(&misc_.kmpo);

  for(unsigned i=0;i<1000;i++)
    shower_(&iqi,&ei,&xi,&yi,&zi,&ui,&vi,&wi,&iri,&wti);

  counters_out_(&ione);
}

void ausgab_(int* iarg)
{
}

void howfar_(void)
{
}

void randomset_(double* d, int* icall)
{
  origrandomset_(d);
}
