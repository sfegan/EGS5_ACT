// BFfield.hpp - Magnetic field as a function of height
// Stephen Fegan - sfegan@llr.in2p3.fr - September 2012
// $Id: BField.hpp 4546 2012-10-01 12:04:16Z sfegan $

#ifndef BFIELD_HPP
#define BFIELD_HPP

class BField
{
public:
  virtual ~BField();
  virtual void getFieldCGS(double h, double& bx, double& by, double& bz) = 0;
};

class ConstantBField: public BField
{
public:
  ConstantBField(double bx, double by, double bz):
    BField(), m_bx(bx), m_by(by), m_bz(bz) { /* nothing to see here */ }
  virtual ~ConstantBField();
  virtual void getFieldCGS(double h, double& bx, double& by, double& bz);
private:
  double m_bx;
  double m_by;
  double m_bz;
};

#endif // defined BFIELD_HPP
