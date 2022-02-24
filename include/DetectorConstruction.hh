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
    // For setting options
    // *******************

    char ldmat; // if 1, LD material is Silicon,
    		// else Germanium

    char emat; 	// if 1, EMMaterial is G4_Cu, 
    		// 2 then Ge_Pyrex_Glass
		// 3 then G4_GLASS_LEAD
		// 4 then G4_POLYTRIFLUOROCHLOROETHYLENE
		// else G4_Al

    char none; 	// if 1, no encapsulatory materials -- only LDs,
   	       	// else, only two LDs
		// overrides em

    char pg;	// if 1, using photon gun

    char em; 	// if 1, using encapsulatory materials 
		// none overrides em; used if none != 1

    char sf;	// if 3, all sides polished
  		// 2 then all sides ground 
		// else, two sides ground, others polished 

    char lmos; 	// if 1, LMO exists
   		// if 2, LMO exists and is a sensitive detector 
		// else, not

    char oc;	// if 1, in optical contact
    		// else, not

    char arc;	// if 1 and oc != 1, use antireflective coating
    		// else, not

    double arct;// thickness of antireflective coating 
    		// if in use (in nms)
    
    char dddm; // if 1, dielectric_dielectric boundaries for EM, LD
    		// if 2, dielectric_metal
		// else no trans, refl defined (using d_d)

    char pLD5; // if 1, pointing at LD5,
    		// turn it into EM
		// (used if none = '1')

    float ldg;	// value of LD gap! (in mm)
    float emg;	// value of EM gap!


    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4int   fNofLayers;     // number of layers
    DetectorMessenger* fDetectorMessenger;

    G4double fLMO_xy;

    // *****
    // World
    // *****

    G4LogicalVolume* worldLV;
    G4MaterialPropertiesTable* fWorldMPT;
    G4Material* WorldMaterial;

    // ***********
    // LMO Crystal
    // ***********

    G4LogicalVolume* fLMO_LV_primary; // LMO_primary
    G4MaterialPropertiesTable* fLMOMPT;
    G4Material* LMOMaterial;
    
    // **************
    // Light Detector
    // **************
    
    G4MaterialPropertiesTable* fLDMPT;
    G4Material* LightDetectorMaterial;

    // **********************
    // Encapsulatory Material 
    // **********************
    
    G4LogicalVolume* EMLV; // EMLV
    G4MaterialPropertiesTable* EMMPT;
    G4Material* EMMaterial;

    // **********************
    // Optical Surfaces
    // **********************

    G4MaterialPropertiesTable* LDk;
    G4MaterialPropertiesTable* EMk;
    //G4MaterialPropertiesTable* dmpP;
    //G4MaterialPropertiesTable* dmgP;
   

    // **********************
    // Antireflective Coating
    // **********************
   
    //G4LogicalVolume* ARCLV; 
    G4MaterialPropertiesTable* ARCMPT;
    G4Material* ARCMaterial;
};

#endif
