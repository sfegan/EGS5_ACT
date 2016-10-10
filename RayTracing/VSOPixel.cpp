//-*-mode:c++; mode:font-lock;-*-

#include <sstream>
#include <cassert>

#include <VSDataConverter.hpp>

#include "VSOPixel.hpp"
#include "VSOTelescope.hpp"

using namespace Physics;
using namespace VERITAS;

VSOPixel::VSOPixel():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos()
{
  // nothing to see here
}

VSOPixel::VSOPixel(const VSOTelescope* T, 
		   unsigned PID, unsigned PHID, bool REM,
		   const Physics::Vec3D& P):
  fTelescope(T), fID(PID), fHexID(PHID), fRemoved(REM), fPos(P)
{
  // nothing to see here
}

VSOPixel::~VSOPixel()
{
  // nothing to see here
}

Vec3D VSOPixel::incomingSkyVectorAtZenith(double plate_scale) const
{
  double x = fPos.x*plate_scale;
  double y = -fPos.z*plate_scale;
  double z = -(fTelescope->focalPlanePosition().y+fPos.y);

  Vec3D p(x,y,z);
  p/=p.Norm();

  return p;
}

void VSOPixel::dumpShort(std::ostream& stream) const
{
  stream
    << "PIXEL "

    << VSDataConverter::toString(fTelescope->id()) << ' '
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fHexID) << ' '
    << VSDataConverter::toString(fRemoved) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << std::endl;
}

void VSOPixel::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fTelescope->id(), "", 30, l) << std::endl
    << FDAV("Pixel ID", fID, "", 30, l) << std::endl
    << FDAV("Pixel Hex ID", fHexID, "", 30, l) << std::endl
    << FDAV("Removed", fRemoved, "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "", 30, l) << std::endl;
}

VSOPixel* VSOPixel::createFromShortDump(std::istream& stream,
					  const VSOTelescope* T)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOPixel* pixel = new VSOPixel;
  pixel->fTelescope=T;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("PIXEL"));

  unsigned TID;

  linestream
    >> TID
    >> pixel->fID
    >> pixel->fHexID
    >> pixel->fRemoved
    >> pixel->fPos.x

    >> pixel->fPos.y
    >> pixel->fPos.z;

  if(!linestream)
    {
      delete pixel;
      return 0;
    }

  return pixel;
}
