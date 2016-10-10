//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOArrayParameters.hpp

  Array parameters class header file.

  AUTOMATICALLY GENERATED BY code_gen_array_ini_hpp.pl DO NOT EDIT!!!

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    11/30/2004
  \version 0.2

*/

#ifndef VSOARRAYPARAMETERS_HPP
#define VSOARRAYPARAMETERS_HPP

#include <iostream>
#include <string>

#include <Vec3D.hpp>

namespace VERITAS
{

  //! Class encapsulates the parameters from which an array is generated.
  class VSOArrayParameters
  {
   public:
    // ************************************************************************
    // Array Parameters
    // ************************************************************************
    double            ArrayLatitude;                      // rad
    double            ArrayLongitude;                     // rad
    double            ArrayAltitude;                      // cm
    unsigned          ArrayNumTelRings;                   // 
    double            ArrayTelSpacing;                    // cm
    bool              ArrayLabelingParity;                // 
    std::string       ScopeMissingList;                   // 
    double            ScopePosXYDisp;                     // cm
    double            ScopePosZDisp;                      // cm
    std::string       ScopePosNorth;                      // cm
    std::string       ScopePosEast;                       // cm
    std::string       ScopePosASL;                        // cm
    double            ScopeDeltaY;                        // rad
    double            ScopeAlphaX;                        // rad
    double            ScopeAlphaY;                        // rad
    double            ScopeElevation;                     // rad
    double            ScopeAzimuth;                       // rad
    double            ScopeTranslationX;                  // cm
    double            ScopeTranslationY;                  // cm
    double            ScopeTranslationZ;                  // cm
    double            ReflectorCurvR;                     // cm
    double            ReflectorApert;                     // cm
    double            ReflectorRot;                       // rad
    unsigned          ReflectorAlignMode;                 // 
    double            ReflectorAlignObjectPlane;          // cm
    double            ReflectorAlignImagePlane;           // cm
    double            ReflectorAlignPtX;                  // cm
    double            ReflectorAlignPtY;                  // cm
    double            ReflectorAlignPtZ;                  // cm
    double            ReflectorAlignPtXZDisp;             // cm
    double            ReflectorAlignPtYDisp;              // cm
    double            ReflectorFPAlignTheta;              // rad
    double            ReflectorFPAlignPhi;                // rad
    double            ReflectorIP;                        // cm
    std::string       ReflectorObscuration;               // 
    double            MirrorFacetSpacing;                 // cm
    double            MirrorFacetSize;                    // cm
    double            MirrorFLength;                      // cm
    double            MirrorFLengthDisp;                  // cm
    double            MirrorSpotSizePhotonFraction;       // 
    double            MirrorSpotSize;                     // cm
    double            MirrorSpotSizeDisp;                 // cm
    double            MirrorDegradingFactor;              // 
    double            MirrorAlignTangentDisp;             // cm
    double            MirrorPosTangentDisp;               // cm
    double            MirrorPosNormalDisp;                // cm
    bool              MirrorLabelingParity;               // 
    std::string       MirrorMissingList;                  // 
    double            CameraDiameter;                     // cm
    double            CameraFPTransX;                     // cm
    double            CameraFPTransY;                     // cm
    double            CameraFPTransZ;                     // cm
    double            CameraFPRotX;                       // rad
    double            CameraFPRotY;                       // rad
    double            CameraFPRotZ;                       // rad
    double            CameraIP;                           // cm
    double            PixelSpacing;                       // cm
    double            PixelConcSurvProb;                  // 
    double            PixelDiameter;                      // cm
    bool              PixelLabelingParity;                // 
    std::string       PixelMissingList;                   // 

    // ************************************************************************
    // Simple Member Functions
    // ************************************************************************
    VSOArrayParameters(const VSOArrayParameters& o);
    VSOArrayParameters(bool use_canonical_values=false);
    virtual ~VSOArrayParameters();
    void dump(std::ostream& stream);
    void reset(bool use_canonical_values=false);

    // ************************************************************************
    // Command Line Member Functions
    // ************************************************************************
    static void zeroCanonicalValues();
    static void resetCanonicalValues();
    static void setCanonicalValuesFromArrayParameters(const VSOArrayParameters& o);

    // ************************************************************************
    // ArrayINI Member Functions
    // ************************************************************************
    bool readFromArrayINIFile(std::istream& stream);
    bool readFromArrayINIFile(const std::string& filename);
    void writeToArrayINIFile(std::ostream& stream) const;
    void writeToArrayINIFile(const std::string& filename) const;

    static const std::string scCollection;

   private:
    static double     sCanonicalArrayLatitude;            // rad
    static double     sCanonicalArrayLongitude;           // rad
    static double     sCanonicalArrayAltitude;            // cm
    static unsigned   sCanonicalArrayNumTelRings;         // 
    static double     sCanonicalArrayTelSpacing;          // cm
    static bool       sCanonicalArrayLabelingParity;      // 
    static std::string sCanonicalScopeMissingList;         // 
    static double     sCanonicalScopePosXYDisp;           // cm
    static double     sCanonicalScopePosZDisp;            // cm
    static std::string sCanonicalScopePosNorth;            // cm
    static std::string sCanonicalScopePosEast;             // cm
    static std::string sCanonicalScopePosASL;              // cm
    static double     sCanonicalScopeDeltaY;              // rad
    static double     sCanonicalScopeAlphaX;              // rad
    static double     sCanonicalScopeAlphaY;              // rad
    static double     sCanonicalScopeElevation;           // rad
    static double     sCanonicalScopeAzimuth;             // rad
    static double     sCanonicalScopeTranslationX;        // cm
    static double     sCanonicalScopeTranslationY;        // cm
    static double     sCanonicalScopeTranslationZ;        // cm
    static double     sCanonicalReflectorCurvR;           // cm
    static double     sCanonicalReflectorApert;           // cm
    static double     sCanonicalReflectorRot;             // rad
    static unsigned   sCanonicalReflectorAlignMode;       // 
    static double     sCanonicalReflectorAlignObjectPlane; // cm
    static double     sCanonicalReflectorAlignImagePlane; // cm
    static double     sCanonicalReflectorAlignPtX;        // cm
    static double     sCanonicalReflectorAlignPtY;        // cm
    static double     sCanonicalReflectorAlignPtZ;        // cm
    static double     sCanonicalReflectorAlignPtXZDisp;   // cm
    static double     sCanonicalReflectorAlignPtYDisp;    // cm
    static double     sCanonicalReflectorFPAlignTheta;    // rad
    static double     sCanonicalReflectorFPAlignPhi;      // rad
    static double     sCanonicalReflectorIP;              // cm
    static std::string sCanonicalReflectorObscuration;     // 
    static double     sCanonicalMirrorFacetSpacing;       // cm
    static double     sCanonicalMirrorFacetSize;          // cm
    static double     sCanonicalMirrorFLength;            // cm
    static double     sCanonicalMirrorFLengthDisp;        // cm
    static double     sCanonicalMirrorSpotSizePhotonFraction; // 
    static double     sCanonicalMirrorSpotSize;           // cm
    static double     sCanonicalMirrorSpotSizeDisp;       // cm
    static double     sCanonicalMirrorDegradingFactor;    // 
    static double     sCanonicalMirrorAlignTangentDisp;   // cm
    static double     sCanonicalMirrorPosTangentDisp;     // cm
    static double     sCanonicalMirrorPosNormalDisp;      // cm
    static bool       sCanonicalMirrorLabelingParity;     // 
    static std::string sCanonicalMirrorMissingList;        // 
    static double     sCanonicalCameraDiameter;           // cm
    static double     sCanonicalCameraFPTransX;           // cm
    static double     sCanonicalCameraFPTransY;           // cm
    static double     sCanonicalCameraFPTransZ;           // cm
    static double     sCanonicalCameraFPRotX;             // rad
    static double     sCanonicalCameraFPRotY;             // rad
    static double     sCanonicalCameraFPRotZ;             // rad
    static double     sCanonicalCameraIP;                 // cm
    static double     sCanonicalPixelSpacing;             // cm
    static double     sCanonicalPixelConcSurvProb;        // 
    static double     sCanonicalPixelDiameter;            // cm
    static bool       sCanonicalPixelLabelingParity;      // 
    static std::string sCanonicalPixelMissingList;         // 
  };
}

#endif // VSOARRAYPARAMETERS_HPP
