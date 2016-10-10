// Atmposphere.cpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: Atmosphere.cpp 4710 2012-11-06 15:00:55Z sfegan $

#include<cmath>
#include<cctype>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<algorithm>

#include"Atmosphere.hpp"

Atmosphere::~Atmosphere()
{
  // nothing to see here
}

void Atmosphere::makeAtmSlices(std::vector<AtmSlice>& atm, unsigned nslice)
{
  makeAtmSlices(atm,nslice,this->topOfAtmosphere(),0);
}

void Atmosphere::makeAtmSlices(std::vector<AtmSlice>& atm, unsigned nslice, 
			       double zmax, double zmin)
{
  atm.resize(nslice);

  double tmax = this->thickness(zmin);
  double tmin = this->thickness(zmax);
  double dt = tmax-tmin;

  double z = zmin;
  for(unsigned islice=0;islice<nslice;islice++)
    {
      AtmSlice& s(atm[islice]);
      double t = dt * double(nslice-islice)/double(nslice);
      s.zb = z;
      s.tb = t;
      t = dt * double(nslice-islice-1)/double(nslice);
      z = this->zForThickness(t+tmin);
      s.zt = z;
      s.tt = t;
      s.rho = (s.tb-s.tt)/(s.zt-s.zb);
    }
}

// ****************************************************************************
// ISOTHERMAL ATMOSPHERE
// ****************************************************************************

IsothermalAtmosphere::
IsothermalAtmosphere(double rho0, double zs, double zmax, double nmo0):
  Atmosphere(), m_ttoa(rho0*zs*std::exp(-zmax/zs)), m_rho0(rho0),
  m_zs(zs), m_zmax(zmax), m_nmo0(nmo0) 
{
  // nothing to see here
}

IsothermalAtmosphere::~IsothermalAtmosphere()
{
  // nothing to see here
}

double IsothermalAtmosphere::rho(double z)
{
  return m_rho0*std::exp(-z/m_zs);
}

double IsothermalAtmosphere::thickness(double z)
{
  return m_rho0*m_zs*std::exp(-z/m_zs)-m_ttoa;
}

double IsothermalAtmosphere::nMinusOne(double z)
{
  return m_nmo0*std::exp(-z/m_zs);
}

double IsothermalAtmosphere::propagationTimeCorrection(double z)
{
  return m_zs*m_nmo0*(1.0-exp(-z/m_zs));
}

double IsothermalAtmosphere::zForThickness(double t)
{
  return -std::log((t+m_ttoa)/(m_rho0*m_zs))*m_zs;
}

double IsothermalAtmosphere::topOfAtmosphere()
{
  return m_zmax;
}

// ****************************************************************************
// LAYERED ATMOSPHERE
// ****************************************************************************

LayeredAtmosphere::LayeredAtmosphere(const std::string& filename):
  Atmosphere(), m_ztoa(), m_ttoa(), m_levels(), m_layers(), m_ilayer()
{
  std::ifstream stream(filename.c_str());
  if(!stream)
    throw std::string("LayeredAtmosphere: could not open: ")+filename;

  std::string line;
  while(std::getline(stream, line))
    {
      unsigned ichar=0;
      while(isspace(line[ichar]))ichar++;
      if(line[ichar] == '#')continue;
      std::istringstream lstream(line);
      Level l;
      lstream >> l.z >> l.rho >> l.t >> l.nmo;
      l.z *= 1e5;
      m_levels.push_back(l);
    }
  initialize();
}

LayeredAtmosphere::LayeredAtmosphere(const std::vector<Level> levels):
  Atmosphere(), m_ztoa(), m_ttoa(), m_levels(levels), m_layers(), m_ilayer()
{
  initialize();
}

template<typename InputIterator>
LayeredAtmosphere::LayeredAtmosphere(InputIterator first, InputIterator last):
  Atmosphere(), 
  m_ztoa(), m_ttoa(), m_levels(first,last), m_layers(), m_ilayer()
{
  initialize();
}

void LayeredAtmosphere::initialize()
{
  if(m_levels.size()<2)
    throw std::string("LayeredAtmosphere: A minimum of 2 levels required.");
  std::sort(m_levels.begin(), m_levels.end(), Level::CmpZAsc());
  m_layers.resize(m_levels.size()-1);

  m_ztoa = m_levels.back().z;
  m_ttoa = m_levels.back().t;
  
  if(m_ttoa <= 0.0)
    {
      // If the thickness at the top of the atmosphere is zero
      // (i.e. the thickness from infinity to that point has been
      // subtracted) then solve for the thickness there by assuming
      // that the scale height between the final three levels is
      // constant

      if(m_levels.size()<3)
	throw std::string("LayeredAtmosphere: A minimum of 3 levels required "
			  "to solve for thickness.");
      
      const double y2 = m_levels[m_levels.size() - 2].t;
      const double y3 = m_levels[m_levels.size() - 3].t;
      const double x1 = m_levels[m_levels.size() - 1].z;
      const double x2 = m_levels[m_levels.size() - 2].z;
      const double x3 = m_levels[m_levels.size() - 3].z;

      double H = (x2-x3)/(std::log(y3)-std::log(y2)); // initial guess
      double num = std::exp(-x3/H)-std::exp(-x1/H);
      double den = std::exp(-x2/H)-std::exp(-x1/H);
      double df = y3/y2 - num/den;
      unsigned niter = 10;
      while(std::fabs(df)/(y3/y2) > 1e-8)
	{
	  // Newton-Ralphson to find value of H giving agreement with data
	  if(niter-- == 0)
	    throw std::string("LayeredAtmosphere: max number of iterations exceeded");
	  double dfdH = 
	    (den*(x3*std::exp(-x3/H)-x1*std::exp(-x1/H))
	     - num*(x2*std::exp(-x2/H)-x1*std::exp(-x1/H)))/(den*den*H*H);
	  H += df/dfdH;
	  num = std::exp(-x3/H)-std::exp(-x1/H);
	  den = std::exp(-x2/H)-std::exp(-x1/H);
	  df = y3/y2 - num/den;
	}

      double t0 = y3/(std::exp(-x3/H)-std::exp(-x1/H));
      m_ttoa = t0*std::exp(-x1/H);

      for(unsigned ilevel=0;ilevel<m_levels.size();ilevel++)
	m_levels[ilevel].t += m_ttoa;
    }

  double ptc = 0;
  for(unsigned ilayer=0;ilayer<m_layers.size();ilayer++)
    {
      Layer& layer(m_layers[ilayer]);
      const Level& t(m_levels[ilayer+1]);
      const Level& b(m_levels[ilayer]);
      layer.zb    = b.z;
      layer.zt    = t.z;
      layer.rhozs = -(t.z-b.z)/(std::log(t.rho)-std::log(b.rho));
      layer.rho0  = b.rho/std::exp(-b.z/layer.rhozs);
      layer.tzs   = -(t.z-b.z)/(std::log(t.t)-std::log(b.t));
      layer.t0    = b.t/std::exp(-b.z/layer.tzs);
      layer.tb    = b.t;
      layer.tt    = t.t;
      layer.nmozs = -(t.z-b.z)/(std::log(t.nmo)-std::log(b.nmo));
      layer.nmo0  = b.nmo/std::exp(-b.z/layer.nmozs);
      ptc += layer.nmozs*layer.nmo0*std::exp(-b.z/layer.nmozs);
      layer.ptc0   = ptc;
      ptc -= layer.nmozs*layer.nmo0*std::exp(-t.z/layer.nmozs);
    }
  m_ilayer = m_layers.end()-1;

  
}

LayeredAtmosphere::~LayeredAtmosphere()
{
  // nothing to see here
}

double LayeredAtmosphere::rho(double z)
{
  findZ(z);
  return m_ilayer->rho0 * std::exp(-z/m_ilayer->rhozs);
}

double LayeredAtmosphere::thickness(double z)
{
  findZ(z);
  return m_ilayer->t0 * std::exp(-z/m_ilayer->tzs) - m_ttoa;
}

double LayeredAtmosphere::nMinusOne(double z)
{
  findZ(z);
  return m_ilayer->nmo0 * std::exp(-z/m_ilayer->nmozs);
}

double LayeredAtmosphere::propagationTimeCorrection(double z)
{
  findZ(z);
  return m_ilayer->ptc0 - 
    m_ilayer->nmozs * m_ilayer->nmo0 * std::exp(-z/m_ilayer->nmozs);
}

double LayeredAtmosphere::zForThickness(double t)
{
  t += m_ttoa;
  std::vector<Layer>::const_iterator ilayer
    = std::lower_bound(m_layers.begin(),m_layers.end(), t, Layer::CmpTDec());
  if(ilayer == m_layers.end())ilayer--;
  return -std::log(t/ilayer->t0)*ilayer->tzs;
}

double LayeredAtmosphere::topOfAtmosphere()
{
  return m_ztoa;
}
