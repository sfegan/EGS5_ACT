// MultiRNG.cpp - Multiple RNGs
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2013
// $Id: MultiRNG.cpp 5698 2013-10-02 13:10:45Z sfegan $

#include "MultiRNG.hpp"

ConstantDeviateSimpleRNG::~ConstantDeviateSimpleRNG()
{
  // nothing to see here
}

double ConstantDeviateSimpleRNG::uniform()
{
  return m_x;
}

PredefinedDeviateSimpleRNG::~PredefinedDeviateSimpleRNG()
{
  // nothing to see here
}

double PredefinedDeviateSimpleRNG::uniform()
{
  double d;
  if(m_ipd!=m_pd.end() && m_session_count==m_ipd->first)d = (m_ipd++)->second;
  else d = m_rng->uniform();
  m_session_count++;
  return d;
}

PredefinedDeviateSimpleRNG::PredefinedDeviateSimpleRNG(SimpleRNG* rng):
  SimpleRNG(), m_rng(rng), m_session_count(), m_total_count(), m_pd(), 
  m_ipd(m_pd.begin())
{
  // nothing to see here  
}

void PredefinedDeviateSimpleRNG::setPD(uint64_t ideviate, double deviate)
{
  m_pd[ideviate]=deviate;
  resetSession();
}

void PredefinedDeviateSimpleRNG::clearPD(uint64_t ideviate)
{
  m_pd.erase(ideviate);
  resetSession();
}

void PredefinedDeviateSimpleRNG::resetSession() 
{
  m_total_count += m_session_count;
  m_session_count=0; 
  m_ipd = m_pd.begin();
}

