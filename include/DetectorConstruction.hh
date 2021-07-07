//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4RunManager.hh"


class G4VPhysicalVolume;
// class G4GlobalMagFieldMessenger;
class DetectorMessenger;

/// Detector construction class to define materials and geometry.
/// The LightDetector is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the LightDetector :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the LightDetector (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of LightDetectorSD type
/// are created and associated with the Absorber and Gap volumes.
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    G4double	GetLMOSize()	{ return LMO_halflength; }
    G4double 	GetLMOXSize() 	{ return fLMO_xy; }
    G4double	GetLMOSep()	{ return SepDist; }
    G4double	GetLDSep()	{ return LDSepDistTop; }

    G4OpticalSurface* GetSurface(void) { return fSurface; }

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void SetWorldMaterial(const G4String&);
    G4Material* GetWorldMaterial() const { return WorldMaterial; }
    void SetLMOMaterial(const G4String&);
    G4Material* GetLMOMaterial() const { return LightDetectorMaterial; }
    void SetLDMaterial(const G4String&);
    G4Material* GetLDMaterial() const { return LMOMaterial; }

    void AddLDMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
    void AddLDMPC(const G4String& prop, G4double v);
    G4MaterialPropertiesTable* GetLDMaterialPropertiesTable()
    {
      return fLDMPT;
    }

    void AddLMOMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
    void AddLMOMPC(const G4String& prop, G4double v);
    G4MaterialPropertiesTable* GetLMOMaterialPropertiesTable()
    {
      return fLMOMPT;
    }

    void AddWorldMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
    void AddWorldMPC(const G4String& prop, G4double v);
    G4MaterialPropertiesTable* GetWorldMaterialPropertiesTable()
    {
      return fWorldMPT;
    }

  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // data members
    //
    // static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger


    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    DetectorMessenger* fDetectorMessenger;
    G4int   fNofLayers;     // number of layers
    G4MaterialPropertiesTable* fLMOMPT;
    G4MaterialPropertiesTable* fLDMPT;
    G4MaterialPropertiesTable* fWorldMPT;
    G4MaterialPropertiesTable* fSurfaceMPT;
    G4Material* WorldMaterial;
    G4Material* LMOMaterial;
    G4Material* LightDetectorMaterial;
    G4LogicalVolume* worldLV;
    G4LogicalVolume* fLMO_LV;
    G4LogicalVolume* LightDetectorLV;


    G4OpticalSurface* fSurface;
    G4double LMO_halflength;
    G4double LD_halfwidth;
    G4double SepDist; // Separation distance between crystals, face-to-face
    G4double SepDistZ; // Separation distance between crystals, face-to-face, z-dir
    G4double SmallSepDist; // Smaller separation distance between crystals, face-to-face
    G4double LDSepDistTop; // Separation distance between crystal & TOP LD face-to-face
    G4double LDSepDistBot; // Separation distance between crystal & TOP LD face-to-face

    G4double fLMO_xy;
    G4double fLD_xy;
    G4double fLD_z;




};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
