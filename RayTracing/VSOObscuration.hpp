//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOObscuration.hpp
  Obscuration class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@llr.in2p3.fr         \n

  \date    07/04/2013
  \version 0.1
  \note
*/

#ifndef VSOOBSCURATION_HPP
#define VSOOBSCURATION_HPP

#include <iostream>
#include <string>
#include <vector>

#include <Vec3D.hpp>
#include <Particle.hpp>

namespace VERITAS
{

  class VSOObscuration
  {
  public:
    virtual ~VSOObscuration();
    virtual bool doesObscure(const Physics::Particle& p_in)
    {
      Physics::Particle p_out;
      return doesObscure(p_in, p_out);
    }
    virtual bool doesObscure(const Physics::Particle& p_in,
			     Physics::Particle& p_out) const = 0;

    virtual std::string dumpToString() const = 0;
    virtual VSOObscuration* clone() const = 0;

    static std::vector<VSOObscuration*> 
    createObsVecFromString(const std::string &str);
    static std::string 
    dumpObsVecToString(const std::vector<VSOObscuration*>& vo);

    static bool tokenize(std::string& str, const std::string& name,
			 std::vector<std::string>& tokens);
  };

  class VSODiskObscuration: public VSOObscuration
  {
  public:
    VSODiskObscuration(const Physics::Vec3D& center,
		       const Physics::Vec3D& normal,
		       double radius, bool incoming_only):
      VSOObscuration(), fX0(center), fN(normal), fR(radius), fD0(),
      fICO(incoming_only)
    {
      fD0 = center*normal;
    }
    virtual ~VSODiskObscuration();
    virtual bool doesObscure(const Physics::Particle& p_in,
			     Physics::Particle& p_out) const;
    virtual std::string dumpToString() const;
    virtual VSOObscuration* clone() const;
    static VSODiskObscuration* createFromString(std::string& str);

  private:
    Physics::Vec3D fX0;
    Physics::Vec3D fN;
    double         fR;
    double         fD0;
    bool           fICO;
  };

  class VSOCylinderObscuration: public VSOObscuration
  {
  public:
    VSOCylinderObscuration(const Physics::Vec3D& x1,
			   const Physics::Vec3D& x2,
			   double radius, bool incoming_only):
      VSOObscuration(), fX1(x1), fX2(x2), fR(radius), fN(), fD1(), fD2(),
      fICO(incoming_only)
    { 
      fN = x2-x1;
      fN /= fN.Norm();
      fD1 = fN*x1;
      fD2 = fN*x2;
      fD  = std::fabs(fD2-fD1);
    }

    virtual ~VSOCylinderObscuration();
    virtual bool doesObscure(const Physics::Particle& p_in,
			     Physics::Particle& p_out) const;
    virtual std::string dumpToString() const;
    virtual VSOObscuration* clone() const;
    static VSOCylinderObscuration* createFromString(std::string& str);

  private:
    Physics::Vec3D fX1;
    Physics::Vec3D fX2;
    double         fR;
    Physics::Vec3D fN;
    double         fD1;
    double         fD2;
    double         fD;
    bool           fICO;
  };
  
}

#endif // defined VSOOBSCURATION
