#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4RunManager.hh"


class G4VPhysicalVolume;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    G4double 	GetLMOXSize() 	{ return fLMO_xy; }
    
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    private:
    
    void DefineMaterials();
    G4VPhysicalVolume* DefineWorld();
    void DefineVolumes();



    // *******************
    // For defining materials
    // *******************

    void DefineLDM(char mat = '1');
    void DefineEMM(char mat = '1');
    void DefineCrysM(char mat = '1');
    void DefineARCM(char mat = '1');
    void DefineWorldM(char mat = '1');

    void DefineMisc(G4Material** CuMaterial, G4Material** AlMaterial);
    void DefineTESMisc(G4Material** FrameMaterial, G4Material** LDMaterial, G4Material** ARCMaterial, G4Material** ReflectorMaterial);

    // *******************
    // For defining physical objects
    // *******************

    G4VPhysicalVolume* DefineCrys(G4double LMOHalfSizes[]);

    void DefineIV(G4double IVHalfSizes[],G4double IVPos[],
G4VPhysicalVolume** IV1,G4VPhysicalVolume** IV2,G4VPhysicalVolume** IV3,G4VPhysicalVolume** IV4,G4VPhysicalVolume** IV5,G4VPhysicalVolume** IV6);

    void DefineLD(G4double LDHalfSizes[],G4double LDPos[],
G4VPhysicalVolume** LD1,G4VPhysicalVolume** LD2,G4VPhysicalVolume** LD3,G4VPhysicalVolume** LD4,G4VPhysicalVolume** LD5,G4VPhysicalVolume** LD6);

    void DefineEM(G4double EMHalfSizes[],G4double EMPos[],
G4VPhysicalVolume** EM1,G4VPhysicalVolume** EM2,G4VPhysicalVolume** EM3,G4VPhysicalVolume** EM4,G4VPhysicalVolume** EM5,G4VPhysicalVolume** EM6);

    void DefineARC(G4double ARCHalfSizes[],G4double ARCPos[],
G4VPhysicalVolume** ARC1,G4VPhysicalVolume** ARC2,G4VPhysicalVolume** ARC3,G4VPhysicalVolume** ARC4,G4VPhysicalVolume** ARC5,G4VPhysicalVolume** ARC6);

	
	
    // *******************
    // Setting
    // *******************

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4int   fNofLayers;     // number of layers
    DetectorMessenger* fDetectorMessenger;

    G4double fLMO_xy;


    // *****
    // TODO
    // *****

    G4LogicalVolume* worldLV;
    
    G4Material* WorldMaterial;
    G4Material* LMOMaterial;
    G4Material* LightDetectorMaterial;
    G4Material* EMMaterial;
    G4Material* ARCMaterial;

    // used in borders' MPTs. :)
    G4MaterialPropertiesTable* LDk;
    G4MaterialPropertiesTable* EMk;
    G4MaterialPropertiesTable* ARCk;
    
};

#endif
