// DetectorEfficiency.hpp - Classes to handle efficiency of detector
// - inculding absoption of Cherenkov light along the path
// Stephen Fegan - sfegan@llr.in2p3.fr - November 2012
// $Id: DetectorEfficiency.hpp 5422 2013-06-26 14:01:03Z sfegan $

#ifndef DETECTOREFFICIENCY_HPP
#define DETECTOREFFICIENCY_HPP

// #define ACT_LIGHT_YIELD_TAYLOR_SERIES_IN_LOG

#include <cmath>

#include <VSDataConverter.hpp>

using namespace VERITAS;

typedef triple<double,double,double> yield_t;

inline yield_t& operator+= (yield_t& a, const yield_t& b)
{
  a.first   += b.first;
  a.second  += b.second;
  a.third   += b.third;
  return a;
}

inline yield_t operator+ (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first+b.first, a.second+b.second, a.third+b.third);
}

inline yield_t operator- (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first-b.first, a.second-b.second, a.third-b.third);
}

inline yield_t operator/ (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first/b.first, a.second/b.second, a.third/b.third);
}

inline yield_t operator* (const yield_t& a, const yield_t& b)
{
  return yield_t(a.first*b.first, a.second*b.second, a.third*b.third);
}

inline yield_t operator* (const yield_t& a, double c)
{
  return yield_t(a.first*c, a.second*c, a.third*c);
}

#include "Interpolation1D.hpp"

#ifdef SWIG
%template(InterpLinear1DDouble) Interpolation1D<double, LinearInterpolator<double> >; 
%template(InterpLinear1DYieldT) Interpolation1D<yield_t, LinearInterpolator<yield_t> >; 
%template(YieldT) VERITAS::triple<double,double,double>;
#endif

class TelescopeEfficiency: public InterpLinear1D
{
public:
  TelescopeEfficiency();
  void scaleEff(const InterpLinear1D& eff);
  void scaleEffFromFile(const std::string& filename, 
			double lambda0_nm=180.0, double dlambda_nm=5.0);
};

#if 0
class ACTILYInterpolator
  : private ExpInterpolator<double>, private LinearInterpolator<double>
{
private:
  typedef ExpInterpolator<double> EXP;
  typedef LinearInterpolator<double> LIN;
public:
  inline yield_t interpolate(double x, 
			     double x0, const yield_t& y0,
			     double x1, const yield_t& y1) const
  {    
    yield_t y;
    y.first  = EXP::interpolate(x, x0, y0.first,  x1, y1.first);
    y.second = LIN::interpolate(x, x0, y0.second, x1, y1.second);
    y.third  = LIN::interpolate(x, x0, y0.third,  x1, y1.third);
    return y;
  }

  inline yield_t integrate(double x0, const yield_t& y0,
			   double x1, const yield_t& y1) const
  {
    yield_t y;
    y.first  = EXP::integrate(x0, y0.first,  x1, y1.first);
    y.second = LIN::integrate(x0, y0.second, x1, y1.second);
    y.third  = LIN::integrate(x0, y0.third,  x1, y1.third);
    return y;
  }
};

#ifdef SWIG
%template(Interp1DACT) Interpolation1D<yield_t, ACTILYInterpolator>;
#endif
#endif

class ACTIntegratedLightYield:
  public Interpolation1D<yield_t, LinearInterpolator<yield_t> >
  //public Interpolation1D<yield_t, ExpInterpolator<yield_t> >
  //public Interpolation1D<yield_t, ACTILYInterpolator>
{
public:
  ACTIntegratedLightYield(double w0=0);
  double yield(double h, double w) const;
  double w0() const { return m_w0; }
private:
  double m_w0;
};

class AtmosphericAbsorption
{
public:
  AtmosphericAbsorption(const std::string& filename, double spacing_km=1.0);
  InterpLinear1D absorptionForAltitude(double h) const;
  ACTIntegratedLightYield integrateYield(double h0, double w0,
					 const TelescopeEfficiency& eff);
private:
  std::vector<double>         m_e_ev;
  std::vector<InterpLinear1D> m_absorption;
};

#endif // not defined DETECTOREFFICIENCY_HPP
