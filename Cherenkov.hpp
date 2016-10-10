// Cherenkov.hpp - Classes to do some simple cherenkov calculations
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: Cherenkov.hpp 5416 2013-06-24 13:46:00Z sfegan $

#ifndef CHERENKOV_HPP
#define CHERENKOV_HPP

#include<cmath>

// Units:
// Energy:    MeV
// Distance:  cm
// Yield:     ph/eV

class Cherenkov
{
public:
  static double c()
  {
    return 2.9979246e+10;
  }

  static double rm()
  {
    return 0.510998927567813;
  }

  static double rmsq()
  {
    return 0.261119903975455;
  }

  static double yield_const()
  {
    return 369.81020849958;
  }

  static double threshold(const double n)
  {
    const double b2 = 1.0/(n*n);
    const double g2 = 1.0/(1-b2);
    return std::sqrt(rmsq()*g2);
  }

  static double sin2Thetac(const double e, const double n)
  {
    const double g2   = e*e/rmsq();         // gamma^2
    const double b2   = 1.0 - 1.0/g2;       // beta^2
    return 1.0 - 1.0/(b2*n*n);
  }
  
  static double yieldDensity(const double sin2_thetac, const double ustep)
  {
    return yield_const()*sin2_thetac*ustep;
  }  
};

#endif
