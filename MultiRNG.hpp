// MultiRNG.hpp - Multiple RNGs
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2013
// $Id: MultiRNG.hpp 5698 2013-10-02 13:10:45Z sfegan $

#ifndef MULTIRNG_HPP
#define MULTIRNG_HPP

#include <stdint.h>
#include <map>

#include "SimpleRNG.hpp"

class ConstantDeviateSimpleRNG: public SimpleRNG
{
public:
  ConstantDeviateSimpleRNG(double x): SimpleRNG(), m_x(x) { }
  virtual ~ConstantDeviateSimpleRNG();
  virtual double uniform();
private:
  double m_x;
};

class PredefinedDeviateSimpleRNG: public SimpleRNG
{
public:
  PredefinedDeviateSimpleRNG(SimpleRNG* rng);
  virtual ~PredefinedDeviateSimpleRNG();
  virtual double uniform();
  void setPD(uint64_t ideviate, double deviate);
  void clearPD(uint64_t ideviate);
  void resetSession();
  uint64_t sessionCount() { return m_session_count; }
  uint64_t totalCount() { return m_total_count+m_session_count; }
private:
  SimpleRNG*                           m_rng;
  uint64_t                             m_session_count;
  uint64_t                             m_total_count;
  std::map<uint64_t, double>           m_pd;
  std::map<uint64_t, double>::iterator m_ipd;
};

#endif
