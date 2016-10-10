#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>

#include "EGS5LayeredDetector.hpp"

typedef EGS5LayeredDetector::Media Media;
typedef EGS5LayeredDetector::Layer Layer;

int main()
{
  try
    {
      std::vector<Media> media;
      media.push_back(Media("AIR AT NTP",10000));

      std::vector<Layer> layers(1);
      layers[0].imedia = 0;
      layers[0].zt     = 10000000;
      
      double emax = 1e8;
      EGS5LayeredDetector det(media, layers, 0, emax);
      EGS5System* egs5 = EGS5System::getEGS5System(&det);
      egs5->initializeEGS5();

      extern EGS5_photin photin_;
      extern EGS5_elecin elecin_;
      extern EGS5_useful useful_;
      extern EGS5_media media_;

      int eiq = -1;
      int piq = 1;
      for(double x=2; x<std::log10(emax); x+=0.1)
	{
	  int medium=0;
	  useful_.medium = medium+1;

	  double eig = std::pow(10.0,x);
	  double gle = x*std::log(10.0);
	  int lgle = int(photin_.ge1[medium]*gle + photin_.ge0[medium]);
	  int iextp = 0;
	  assert(eig >= 0.15);
	  double gmfpr0 = 
	    photin_.gmfp1[medium][lgle+iextp-1]*gle +
	    photin_.gmfp0[medium][lgle+iextp-1];
	  double gbr1 = 
	    photin_.gbr11[medium][lgle+iextp-1]*gle +
	    photin_.gbr10[medium][lgle+iextp-1];
	  double gbr2 = 
	    photin_.gbr21[medium][lgle+iextp-1]*gle + 
	    photin_.gbr20[medium][lgle+iextp-1];
	  
	  double eie = eig;
          double eke = eie - useful_.rm;
          double elke = std::log(eke);
	  int lelke = int(elecin_.eke1[medium]*elke + elecin_.eke0[medium]);

	  double esig0 = 0;
          hardx_(&eiq,&eke,&lelke,&elke,&esig0);

	  double psig0 = 0;
          hardx_(&piq,&eke,&lelke,&elke,&psig0);

	  std::cout << x << ' ' << eig << ' ' << gmfpr0 << ' '
		    << gbr1 << ' ' << gbr2-gbr1 << ' ' << 1-gbr2 << ' '
		    << 1.0/esig0 << ' ' << 1.0/psig0 << '\n';
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
}
