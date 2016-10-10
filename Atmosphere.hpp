// Atmposphere.hpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: Atmosphere.hpp 5416 2013-06-24 13:46:00Z sfegan $

#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

#include<vector>
#include<cassert>

// Units:
// Height, z:    cm
// Density, rho: g/cm^3
// Thickness, t: g/cm^2

class AtmSlice
{
public:
  AtmSlice(double _v=0): zb(_v), zt(_v), tb(_v), tt(_v), rho(_v), 
			 n_minus_one(_v) { }
  double zb;
  double zt;
  double tb;
  double tt;
  double rho;
  double n_minus_one;

  bool operator< (const AtmSlice& o) const { return zt < o.zt; }

  class CmpZAsc
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b) 
    { return a.zt < b.zt; } };
  
  class CmpTDec
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b) 
    { return b.tb < a.tb; } };
};

class Atmosphere
{
public:
  virtual ~Atmosphere();
  virtual double rho(double z) = 0;
  virtual double thickness(double z) = 0;
  virtual double nMinusOne(double z) = 0;
  virtual double propagationTimeCorrection(double z) = 0;
  virtual double zForThickness(double t) = 0;
  virtual double topOfAtmosphere() = 0;
  void makeAtmSlices(std::vector<AtmSlice>& atm, unsigned nslice, 
		     double zmax, double zmin);
  void makeAtmSlices(std::vector<AtmSlice>& atm, unsigned nslice);
};

class IsothermalAtmosphere: public Atmosphere
{
public:
  IsothermalAtmosphere(double rho0=1.2e-3, double zs = 8.5e5, 
		       double zmax = 1.2e7, double nmo0=2.75e-4);
  virtual ~IsothermalAtmosphere();
  virtual double rho(double z);
  virtual double thickness(double z);
  virtual double nMinusOne(double z);
  virtual double propagationTimeCorrection(double z);
  virtual double zForThickness(double t);
  virtual double topOfAtmosphere();
private:
  double m_ttoa;
  double m_rho0;
  double m_zs;
  double m_zmax;
  double m_nmo0;
};  

class LayeredAtmosphere: public Atmosphere
{
public:
  struct Level
  {
    double z;
    double rho;
    double t;
    double nmo;

    class CmpZAsc { public: bool operator() (const Level& a, const Level& b) 
      { return a.z < b.z; } };

    class CmpTDec { public: bool operator() (const Level& a, const Level& b) 
      { return b.t < a.t; } };
  };

  LayeredAtmosphere(const std::string& filename);
  LayeredAtmosphere(const std::vector<Level> levels);
  template<typename InputIterator>
  LayeredAtmosphere(InputIterator first, InputIterator last);

  virtual ~LayeredAtmosphere();
  virtual double rho(double z);
  virtual double thickness(double z);
  virtual double nMinusOne(double z);
  virtual double propagationTimeCorrection(double z);
  virtual double zForThickness(double t);
  virtual double topOfAtmosphere();

  const std::vector<Level>& getLevels() const { return m_levels; }
private:
  struct Layer
  {
    Layer(double _v = 0): 
      zb(_v), zt(_v), rho0(_v), rhozs(_v), t0(_v), tzs(_v), tb(_v), tt(_v),
      nmo0(_v), nmozs(_v), ptc0(_v) { }
    double zb;
    double zt;
    double rho0;
    double rhozs;
    double t0;
    double tzs;
    double tb;
    double tt;
    double nmo0;
    double nmozs;
    double ptc0;
    bool operator< (const Layer& o) const { return zt<o.zt; }

    class CmpTDec { public: bool operator() (const Layer& a, const Layer& b) 
      { return b.tt<a.tt; } };
  };

  void initialize();
  inline void findZ(double z) const
  {
    //    std::cout << "A: " << z << ' ' << m_ilayer-m_layers.begin() << ' ' << m_ilayer->zb << ' ' << m_ilayer->zt << '\n';

    if(z<=m_ilayer->zt)
      {
	if(z<=m_ilayer->zb && m_ilayer!=m_layers.begin())
	  {
	    m_ilayer--;
	    if(z<=m_ilayer->zb && m_ilayer!=m_layers.begin())
	      {
		m_ilayer = std::lower_bound(m_layers.begin(), m_ilayer, z);
		if(z<=m_ilayer->zb)m_ilayer = m_layers.begin();
	      }
	  }
      }
    else
      {
	if(++m_ilayer == m_layers.end())m_ilayer--;
	else if(m_ilayer->zt < z)
	  {
	    ++m_ilayer;
	    m_ilayer = std::lower_bound(m_ilayer,m_layers.end(), z);
	    if(m_ilayer == m_layers.end())m_ilayer--;
	  }
      }
    
    //    std::cout << "B: " << z << ' ' << m_ilayer-m_layers.begin() << ' ' << m_ilayer->zb << ' ' << m_ilayer->zt << '\n';

    assert(m_ilayer<m_layers.end());
    assert(z<=m_ilayer->zt || m_ilayer==m_layers.end()-1);
    assert(m_ilayer->zb<z || m_ilayer==m_layers.begin());
  }

  double m_ztoa;
  double m_ttoa;
  std::vector<Level> m_levels;
  std::vector<Layer> m_layers;
  mutable std::vector<Layer>::const_iterator m_ilayer;
};

#endif
