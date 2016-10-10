// EGS5LayeredDetector.cpp - Detector class for EGS5 system that implements
// - layered detector geometry and magnetic field (as in ucbend.f).
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: EGS5LayeredDetector.cpp 4954 2013-01-18 08:44:50Z sfegan $

#include <cassert>
#include <algorithm>
#include <limits>
#include <iostream> 

#include "EGS5LayeredDetector.hpp"

EGS5LayeredDetector::EGS5LayeredDetector(): 
  EGS5UserInterface(), m_media(), m_layers(), m_zbot(), m_emax(), 
  m_rmsq(), m_n_howfar_entry(), m_n_howfar_change_region()
{
  // nothing to see here
}

void EGS5LayeredDetector::
initialize(const std::vector<Media>& media, const std::vector<Layer>& layers, 
	   double zbot, double emax)
{
  m_media = media;
  m_layers.clear();
  m_layers.insert(m_layers.begin(), layers.begin(), layers.end());
  m_zbot = zbot;
  m_emax = emax;
  extern EGS5_useful useful_;
  std::sort(m_layers.begin(), m_layers.end());
  if(m_layers.front().zt < m_zbot)
    throw std::string("EGS5LayeredDetector: zbot must be smaller than \"top\" "
		      "of first layer");
  m_layers.front().zb = m_zbot;  
  for(unsigned ilayer=1;ilayer<m_layers.size();ilayer++)
    m_layers[ilayer].zb = m_layers[ilayer-1].zt;
  m_rmsq = useful_.rm*useful_.rm;
}

EGS5LayeredDetector::
EGS5LayeredDetector(const std::vector<Media>& media, 
		    const std::vector<Layer>& layers, 
		    double zbot, double emax):
  EGS5UserInterface(), m_media(), m_layers(), m_zbot(), m_emax(), 
  m_rmsq(), m_n_howfar_entry(), m_n_howfar_change_region()
{
  initialize(media, layers, zbot, emax);
}

EGS5LayeredDetector::~EGS5LayeredDetector()
{
  // nothing to see here
}

double EGS5LayeredDetector::getEMax()
{
  return m_emax;
}

void EGS5LayeredDetector::getMedia(std::vector<Media>& media)
{
  media = m_media;
}

void EGS5LayeredDetector::getRegions(std::vector<Region>& regions)
{
  regions.resize(m_layers.size() + 2);
  regions.front().imedia = regions.back().imedia = -1;
  regions.front().rho    = regions.back().rho    = 0;  
  for(unsigned ilayer=0;ilayer<m_layers.size();ilayer++)
    {
      regions[ilayer+1].imedia = m_layers[ilayer].imedia;
      regions[ilayer+1].rho    = m_layers[ilayer].rho;
    }
}

void EGS5LayeredDetector::ausgab(int& iarg)
{
  // nothing to see here
}

void EGS5LayeredDetector::howfar(void)
{
  extern EGS5_stack stack_;
  extern EGS5_epcont epcont_;

  m_n_howfar_entry++;

  int ip = stack_.np-1;
  int ianow = stack_.ir[ip]-2;

  if(ianow == -1)
    {
      epcont_.idisc=1;
      return;
    }

  double z(stack_.z[ip]);
  double w(stack_.w[ip]);
  int iq(stack_.iq[ip]);

  //std::cout << "A: " << ianow << ' ' << z << ' ' << w << ' ' << iq << ' ' << stack_.e[ip] << '\n';

  int ianxt = ianow;
  double deltaz = std::numeric_limits<double>::infinity();

  if(ianow == int(m_layers.size()))
    {
      if(w > 0)
	{
	  epcont_.idisc=1;
	  return;
	}
      //stack_.dnear[ip] = z-m_layers.back().zt;
      deltaz = (m_layers.back().zt-z)/w;
      ianxt = m_layers.size()-1;
    }
#if 0
  else if(ianow == -1)
    {
      if(w < 0)
	{
	  epcont_.idisc=1;
	  return;
	}
      //stack_.dnear[ip] = m_layers.front().zb-z;
      deltaz = (m_layers.front().zb-z)/w;
      ianxt = 0;
    }
#endif
  else if((iq==0) || (!m_layers[ianow].bfield))
    {
      const InternalLayer& layer(m_layers[ianow]);
      //stack_.dnear[ip] = std::min(layer.zt-z, z-layer.zb);
      if(w<0)deltaz = (layer.zb-z)/w, ianxt=ianow-1;
      else if(w>0)deltaz = (layer.zt-z)/w, ianxt=ianow+1;
    }
  else
    {
      // Charged particle in magnetic field.. bend the trajectory as
      // in the ucbend example - this is probably not a great way to
      // handle this, from a numerical point of view. It would be
      // preferable to propagate the particle with the velocity at the
      // mid-point of the step, whereas here it is done with the
      // velocity at the end-point.
      const InternalLayer& layer(m_layers[ianow]);
      //stack_.dnear[ip] = 0;
      double e(stack_.e[ip]);
      double ux(stack_.u[ip]);
      double uy(stack_.v[ip]);
      double uz(stack_.w[ip]);
      double bx(layer.ubx);
      double by(layer.uby);
      double bz(layer.ubz);
      double psq = std::max(e*e - m_rmsq, 1e-9);
      double p = std::sqrt(psq);
      double rcurv = E5S_RCURV_CM_PER_MEV_OVER_GAUSS*p/layer.btotal;
      double alpha = epcont_.ustep/rcurv;
#warning Should place limits on alpha and ustep here
      if(iq<0)rcurv = -rcurv, alpha = -alpha;
      double ub = bx*ux + by*uy + bz*uz;
      double dzmax = std::min(z-layer.zb, layer.zt-z);
      if(dzmax < epcont_.ustep)
	{
	  double ustep = std::numeric_limits<double>::infinity();
	  for(unsigned ibdy=0;ibdy<2;ibdy++)
	    {
	      double dz = ibdy?(layer.zt-z):(layer.zb-z);
	      double a = (ux*by-uy*bx)*rcurv-(1-ub*ub)*dz;
	      double b = uz*rcurv;
	      double c = -dz;
	      double disc2 = b*b-4*a*c;
	      if(disc2>=0)
		{
		  double q = -0.5*(b+copysign(1.0,b)*std::sqrt(disc2));
		  double us1 = q/a*rcurv;
		  double us2 = c/q*rcurv;
		  if(us1>0 && us1<ustep)ustep=us1;
		  if(us2>0 && us2<ustep)ustep=us2;
		}
	    }
	  if(ustep < epcont_.ustep)alpha = ustep/rcurv;
	}

      double renorm = 1.0/std::sqrt(1.0+alpha*alpha*(1.0-ub*ub));
      stack_.u[ip]     = (ux + alpha*(uy*bz - uz*by))*renorm;
      stack_.v[ip]     = (uy + alpha*(uz*bx - ux*bz))*renorm;
      stack_.w[ip] = w = (uz + alpha*(ux*by - uy*bx))*renorm;
      if(w<0)deltaz = (layer.zb-z)/w, ianxt=ianow-1;
      else if(w>0)deltaz = (layer.zt-z)/w, ianxt=ianow+1;
    }
  
#if 0
  std::cout << iq << ' ' << m_layers[ianow].bfield << ' ' << z << ' ' 
	    << m_layers[ianow].zb << ' ' << m_layers[ianow].zt << ' '
	    << ianow << ' ' << ianxt << ' ' << deltaz << '\n';
#endif

  if(deltaz<=epcont_.ustep)
    {
      int irnxt = ianxt+2;
      epcont_.ustep = deltaz;
      epcont_.irnew = irnxt;
      m_n_howfar_change_region++;
    }
}

void EGS5LayeredDetector::getLayers(std::vector<Layer>& layers)
{
  layers.resize(m_layers.size());
  for(unsigned ilayer=0;ilayer<layers.size();ilayer++)
    layers[ilayer] = m_layers[ilayer];
}
