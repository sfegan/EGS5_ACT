#include "BField.hpp"

BField::~BField()
{
  // nothing to see here
}

ConstantBField::~ConstantBField()
{
  // nothing to see here
}

void ConstantBField::getFieldCGS(double h, double& bx, double& by, double& bz)
{
  bx = m_bx;
  by = m_by;
  bz = m_bz;
}
