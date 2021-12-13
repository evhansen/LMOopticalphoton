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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "LightDetectorSD.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"

#include "G4PVPlacement.hh"
// #include "G4PVReplica.hh"
// #include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"



#include <stdio.h>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
// G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;



// Aux Functions (darn external macros are impossible to figure out)
char readExtFile(const char* filename);






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fDetectorMessenger(nullptr),
   fNofLayers(-1)
{
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fLMOMPT;
  delete fLDMPT;
  delete fWorldMPT;
  delete fSurfaceMPT;
  delete fSurface;
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();


  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

    // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density, iz;
  G4String name= "";
  G4String symbol = "";
  G4int natoms;

  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density
  G4Element* elLi = new G4Element(name="Lithium", symbol="Li", iz=3., a=6.94*g/mole);
  G4Element* elMo = new G4Element(name="Molybdenum",symbol="Mo", iz=42, a=95.96*g/mole);
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=16.00*g/mole);

  // define LMO
  density = 3.073 * g/cm3; // 280 g, 45mm cube
  G4Material* Li2MoO4 = new G4Material (name="Li2MoO4", density, 3);
  Li2MoO4->AddElement(elLi, natoms=2);
  Li2MoO4->AddElement(elMo, natoms=1);
  Li2MoO4->AddElement(elO, natoms=4);

  // Vacuum
  new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
  // define LD
  density = 5.323*g/cm3;
  a = 72.61*g/mole;
  z = 32.0;
  G4Material* Germanium = new G4Material(name="Germanium"   , z , a, density );

  density = 2.33*g/cm3;
  a = 28.0855*g/mole;
  z = 14.0;
  G4Material* Silicon = new G4Material(name="Silicon"   , z , a, density );

  /*
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(1);

  //density = 2.7*g/cm3; // note at 1 atm
  //a = 26.982*g/mole;
  //z = 13.0;
  G4Material* Aluminium = new G4Material(name="Aluminium",z,a,density);
  
  //G4bool isotope = false;
  //G4Material* Aluminium = manager->FindOrBuildMaterial("G4_Al",isotope);
  // Am going to use the NIST DB
  //  Can check available elements and material via 
  //  /material/nist/printElements,listMaterials
  */

  // MATERIAL PROPERTIES TABLES
  fLMOMPT    = new G4MaterialPropertiesTable();
  fLDMPT    = new G4MaterialPropertiesTable();
  fWorldMPT   = new G4MaterialPropertiesTable();
  fSurfaceMPT = new G4MaterialPropertiesTable();
 
  G4int numentries = 3;
  G4double abslength[3] = {100*mm, 100*mm, 100*mm};
  G4double rindex[3] = {1.44, 1.44, 1.44};

  G4double LDabslength[3] = {0.1*mm, 0.1*mm, 0.1*mm};
  G4double LDrindex[3] = {1.97, 1.97, 1.97};

  // G4double worldabslength[3] = {100*mm, 100*mm, 100*mm};
  G4double energies[3] = {2.0*eV, 3.0*eV, 4.0*eV};

  fLMOMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
  fLMOMPT->AddProperty("RINDEX", energies, rindex, numentries);

  fLDMPT->AddProperty("RINDEX", energies, LDrindex, numentries);
  fLDMPT->AddProperty("ABSLENGTH", energies, LDabslength, numentries);
 
  // *
  // TODO: figure out if we actually need to init materials here? 
  //  Seems like probably not
  // *

  /*
  // we want to keep the properties from the NIST DB
  //G4MaterialPropertiesTable* fARMPT = Aluminium->GetMaterialPropertiesTable();
  G4MaterialPropertiesTable* fARMPT = new G4MaterialPropertiesTable();
  
  G4double AlRR[3] = {1.366,0.52135,0.28003};
  G4double AlIR[3] = {7.4052,5.0008,3.7081};
  //G4double AlAb[3] = {};
  //G4double Ref[3] = {0.98,0.98,0.98};
  // 1 eV ~ 1.2398\times 10^(-6) m ~ 1 (micro)m
  // 2 eV ~ 6.1992\times 10^(-7) m ~ 0.06 (micro)m
  // 3 eV ~ 4.1328\times 10^(-7) m ~ 0.04 (micro)m
  // 4 ev ~ 3.0996\times 10^(-7) m ~ 0.03 (micro)m

  //fARMPT->AddProperty("REFLECTIVITY",energies,Ref,numentries);
  
  //fARMPT->AddProperty("ABSLENGTH",energies,AlAb,numentries);
  fARMPT->AddProperty("REALRINDEX",energies,AlRR,numentries);
  fARMPT->AddProperty("IMAGINARYRINDEX",energies,AlIR,numentries);
  // add in what we need

  Aluminium->SetMaterialPropertiesTable(fARMPT);
  // add the table with changes 
  */

  // fWorldMPT->AddProperty("ABSLENGTH",energies, worldabslength, numentries);

  fSurface = new G4OpticalSurface("Surface");
  fSurface->SetType(dielectric_dielectric);
  fSurface->SetFinish(ground);
  fSurface->SetModel(unified);
  fSurface->SetMaterialPropertiesTable(fSurfaceMPT);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 1;
  // G4double absoThickness = 10.*mm;
  // G4double gapThickness =  5.*mm;
  // G4double LightDetectorSizeXY  = 10.*cm;

  // auto layerThickness = absoThickness + gapThickness;
  // auto LightDetectorThickness = fNofLayers * layerThickness;
  // auto worldSizeXY = 1.2 * LightDetectorSizeXY;
  // auto worldSizeZ  = 1.2 * LightDetectorThickness;


  G4double worldSizeXY = 10. * m;

  // LMO & LD dimensions
  LMO_halflength = 0.5 * 45. * mm;
  LD_halfwidth = 0.5 * 0.5 * mm;

  fLMO_xy = LMO_halflength;
  fLD_xy = LMO_halflength;
  fLD_z = LD_halfwidth;

  SepDist = 15. * mm; // TOTAL DISTANCE, Face to Face
  SmallSepDist =  8 * mm;
  SepDistZ = 3. * mm;
  LDSepDistTop = 2. * mm;
  LDSepDistBot = 0.5 * mm;

  G4double fLMO_C2C_inTower = 0.5*SmallSepDist + LMO_halflength*2.;
  G4double fLMO_C2C_diffTower = 0.5*SepDist + LMO_halflength*2.;
  G4double fLD_C2C_Top = LDSepDistTop + LMO_halflength + LD_halfwidth;
  G4double fLD_C2C_Bot = LDSepDistBot + LMO_halflength + LD_halfwidth;

  G4double fLMO_UpperRowZ = fLD_C2C_Top + fLD_C2C_Bot;
  G4double fLD_UpperRowZ = fLMO_UpperRowZ + fLD_C2C_Top;
  G4double fLD_LowerRowZ = fLMO_UpperRowZ + fLD_C2C_Bot;

  // Get materials
  // WorldMaterial = G4Material::GetMaterial("Vacuum");
  // LightDetectorMaterial = G4Material::GetMaterial("Silicon");
  // LMOMaterial = G4Material::GetMaterial("Li2MoO4");

  WorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Vacuum");
  LightDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Germanium");
  LMOMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Li2MoO4");

  WorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
  LightDetectorMaterial->SetMaterialPropertiesTable(fLDMPT);
  LMOMaterial->SetMaterialPropertiesTable(fLMOMPT);

  // *
  // Since modifying the MPT of G4_Al causes a core dump, just make a new MPT
  //  and assign such.
  // * 
  G4bool isotope = false;
  G4Material* ARMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al",isotope);

  G4MaterialPropertiesTable* fARMPT = new G4MaterialPropertiesTable();
  G4double energies[3] = {2*eV,3*eV,4*eV};
  G4double AlRR[3] = {1.366,0.52135,0.28003};
  G4double AlIR[3] = {7.4052,5.0008,3.7081};
  G4int numentries = 3;
  // 1 eV ~ 1.2398\times 10^(-6) m ~ 1 (micro)m
  // 2 eV ~ 6.1992\times 10^(-7) m ~ 0.06 (micro)m
  // 3 eV ~ 4.1328\times 10^(-7) m ~ 0.04 (micro)m
  // 4 ev ~ 3.0996\times 10^(-7) m ~ 0.03 (micro)m
  fARMPT->AddProperty("REALRINDEX",energies,AlRR,numentries);
  fARMPT->AddProperty("IMAGINARYRINDEX",energies,AlIR,numentries);
  ARMaterial->SetMaterialPropertiesTable(fARMPT);
  

  /* MPT modification that doesn't work
  G4MaterialPropertiesTable* fARMPT = ARMaterial->GetMaterialPropertiesTable();
  G4double energies[3] = {2*eV,3*eV,4*eV};
  G4double AlRR[3] = {1.366,0.52135,0.28003};
  G4double AlIR[3] = {7.4052,5.0008,3.7081};
  G4int numentries = 3;
  // 1 eV ~ 1.2398\times 10^(-6) m ~ 1 (micro)m
  // 2 eV ~ 6.1992\times 10^(-7) m ~ 0.06 (micro)m
  // 3 eV ~ 4.1328\times 10^(-7) m ~ 0.04 (micro)m
  // 4 ev ~ 3.0996\times 10^(-7) m ~ 0.03 (micro)m
  fARMPT->AddProperty("REALRINDEX",energies,AlRR,numentries);
  fARMPT->AddProperty("IMAGINARYRINDEX",energies,AlIR,numentries);
  ARMaterial->SetMaterialPropertiesTable(fARMPT);
  */

  if ( ! WorldMaterial || ! LightDetectorMaterial || ! LMOMaterial || !ARMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeXY/2); // its size

  worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 WorldMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 900,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // LMO
  //
  auto LMO_box
    = new G4Box("LMO",          // its name
                fLMO_xy, fLMO_xy, fLMO_xy);    //its size

  fLMO_LV_primary
    = new G4LogicalVolume(
                LMO_box,         // its solid
                LMOMaterial,    // its material
  		"LMO_primary");          // its name

  //fLMO_LV_primary->SetVisAttributes(G4Colour::Yellow());
  G4VisAttributes* LMOVA = new G4VisAttributes();
  LMOVA->SetColor(1.0,1.0,0.0,0.2);
  LMOVA->SetForceSolid(true);
  fLMO_LV_primary->SetVisAttributes(LMOVA);

  fLMO_LV
    = new G4LogicalVolume(
                LMO_box,         // its solid
                LMOMaterial,    // its material
                "LMO");          // its name

  //fLMO_LV->SetVisAttributes(G4Colour::Red());

  G4VPhysicalVolume* fLMO1 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(),  // place off center eventually
                fLMO_LV_primary,          // its logical volume
                "LMO1",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                101,                // copy number
                fCheckOverlaps);  // checking overlaps

  // *
  // Setting up the necessary rotations for surrounding physical volumes (to LMO)
  // *
  G4RotationMatrix* t1 = new G4RotationMatrix();
  t1->rotateY(90.*deg);
  t1->rotateX(90.*deg);

  G4RotationMatrix* p1 = new G4RotationMatrix();
  p1->rotateY(90.*deg);
  p1->rotateZ(90.*deg);

  // *
  // Interim Vacuum (between LMO and LD)
  // *
  G4double IVT = 2.5*mm; // Interim Vacuum Thickness

  // Overlap
  G4double LMOOL = 0.5*mm;
  G4double LDOL = 2.25*mm; 

  IVT = 1.25*mm-0.125*mm-0.195*mm-0.055*mm;

  // Staggering to find accurate positions by brute force
  G4double IVTS = 0.75*mm+0.125*mm+0.25*mm; // (IVTS \in [0,IVT] for (IVT \in [0,2])

  auto InterimVacuum
    = new G4Box("InterimVacuum",     // its name
                 fLD_xy, fLD_xy, IVT); // its size

  G4LogicalVolume* IVLV
    = new G4LogicalVolume(
                 InterimVacuum,     // its solid
                 WorldMaterial,  // its material
                 "InterimVacuumLV");   // its name

  //G4VisAttributes* IVVA = new G4VisAttributes();
  //IVVA->SetColor(0.0,0.0,0.0,5.0);
  //IVVA->SetForceSolid(true);
  //IVLV->SetVisAttributes(IVVA);

  G4VPhysicalVolume* IV1 = new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,(fLMO_C2C_inTower)*(0.5)-IVTS),  // at (0,0,0) 
	       IVLV,          // its logical volume
               "IV1T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               20,                // copy number
               fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* IV2 = new G4PVPlacement(
               0,                // no rotation 
               G4ThreeVector(0,0,(-fLMO_C2C_inTower)*(0.5)+IVTS),  // at (0,0,0)
	       IVLV,          // its logical volume
               "IV1B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               21,                // copy number
               fCheckOverlaps);  // checking overlaps


  G4VPhysicalVolume* IV3 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(-fLMO_C2C_inTower)*(0.5)+IVTS,0),  // at (0,0,0)
              IVLV,          // its logical volume
              "IV2T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              22,                // copy number
              fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* IV4 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(fLMO_C2C_inTower)*(0.5)-IVTS,0),  // at (0,0,0)
              IVLV,          // its logical volume
              "IV2B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              23,                // copy number
              fCheckOverlaps);  // checking overlaps
 
  G4VPhysicalVolume* IV5 = new G4PVPlacement(
               p1,                // no rotation
               G4ThreeVector((-fLMO_C2C_inTower)*(0.5)+IVTS,0,0),  // at (0,0,0)
               IVLV,          // its logical volume
               "IV3T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               24,                // copy number
               fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* IV6 = new G4PVPlacement(
               p1,                // no rotation
               G4ThreeVector((fLMO_C2C_inTower)*(0.5)-IVTS,0,0),  // at (0,0,0)
               IVLV,          // its logical volume
               "IV3B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               25,                // copy number
               fCheckOverlaps);  // checking overlaps

  // *
  // TODO: logical volumes for all individual light detector and
  //  aluminium reflector instances so that they can colour code
  //  their copy numbers?
  // *

  // *
  // LightDetector
  // *
  
  //G4double excess = 1.5*mm;
 
  auto LightDetector
    = new G4Box("LightDetector",     // its name
                 fLD_xy, fLD_xy, fLD_z); // its size

  LightDetectorLV
    = new G4LogicalVolume(
                 LightDetector,     // its solid
                 LightDetectorMaterial,  // its material
                 "LightDetectorLV");   // its name

  //G4VisAttributes* LDVA = new G4VisAttributes(G4Colour::Blue());
  //LDVA->SetForceLineSegmentsPerCircle(4);
  
  //LightDetectorLV->SetVisAttributes(G4Colour::Blue());
  
  //G4VisAttributes* LDVA = new G4VisAttributes();
  //LDVA->SetColor(0.0,0.0,1.0,1.0);
  //LDVA->SetForceLineSegmentsPerCircle(4);
  //LightDetectorLV->SetVisAttributes(LDVA);

  /*
  G4VPhysicalVolume* LD1 = new G4PVPlacement(
               0,                // no rotation
               //G4ThreeVector(0,0,fLD_C2C_Top),  // at (0,0,0)
               G4ThreeVector(0,0,(fLMO_C2C_inTower)*(0.5)),  // at (0,0,0) 
	       LightDetectorLV,          // its logical volume
               "LD1T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               0,                // copy number
               fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* LD2 = new G4PVPlacement(
               0,                // no rotation 
               //G4ThreeVector(0,0,-fLD_C2C_Bot-excess),  // at (0,0,0)
               G4ThreeVector(0,0,(-fLMO_C2C_inTower)*(0.5)),  // at (0,0,0)
	       LightDetectorLV,          // its logical volume
               "LD1B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               1,                // copy number
               fCheckOverlaps);  // checking overlaps


  G4VPhysicalVolume* LD3 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(-fLMO_C2C_inTower)*(0.5),0),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD2T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              2,                // copy number
              fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* LD4 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(fLMO_C2C_inTower)*(0.5),0),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD2B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              3,                // copy number
              fCheckOverlaps);  // checking overlaps
 */
  
  G4VPhysicalVolume* LD5 = new G4PVPlacement(
               p1,                // no rotation
               G4ThreeVector((-fLMO_C2C_inTower)*(0.5),0,0),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD3T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               0,                // copy number
               fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* LD6 = new G4PVPlacement(
               p1,                // no rotation
               G4ThreeVector((fLMO_C2C_inTower)*(0.5),0,0),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD3B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               1,                // copy number
               fCheckOverlaps);  // checking overlaps


  // *
  // Aluminium Reflectors
  // *
  auto AR
    = new G4Box("AluminiumReflectors",     // its name
                 fLD_xy, fLD_xy, fLD_z); // its size

  G4LogicalVolume* ARLV
    = new G4LogicalVolume(
                 AR,     // its solid
                 ARMaterial,  // its material
                 "ARLV");   // its name

  //G4int segs = 4;
  G4VisAttributes* ARVA = new G4VisAttributes();
  ARVA->SetColor(0.0,1.0,0.0,1.0);
  ARVA->SetForceLineSegmentsPerCircle(4);
  ARLV->SetVisAttributes(ARVA);


  G4VPhysicalVolume* AR1 = new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,(fLMO_C2C_inTower)*(0.5)),  // at (0,0,0) 
	       ARLV,          // its logical volume
               "AR1T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               3,                // copy number
               fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* AR2 = new G4PVPlacement(
               0,                // no rotation 
               G4ThreeVector(0,0,(-fLMO_C2C_inTower)*(0.5)),  // at (0,0,0)
	       ARLV,          // its logical volume
               "AR1B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               4,                // copy number
               fCheckOverlaps);  // checking overlaps


  G4VPhysicalVolume* AR3 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(-fLMO_C2C_inTower)*(0.5),0),  // at (0,0,0)
              ARLV,          // its logical volume
              "AR2T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              5,                // copy number
              fCheckOverlaps);  // checking overlaps

  G4VPhysicalVolume* AR4 = new G4PVPlacement(
              t1,                // 
              G4ThreeVector(0,(fLMO_C2C_inTower)*(0.5),0),  // at (0,0,0)
              ARLV,          // its logical volume
              "AR2B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              6,                // copy number
              fCheckOverlaps);  // checking overlaps 


 
  //
  // Interfaces
  //


  // *
  // Setting up surfaces and surface types for interfaces.
  // *  
  G4LogicalBorderSurface* surface1 = new G4LogicalBorderSurface("Surface1", fLMO1, worldPV, fSurface);
  dynamic_cast<G4OpticalSurface*>(surface1->GetSurface(fLMO1, worldPV)->GetSurfaceProperty());

  G4OpticalSurface* t2 = new G4OpticalSurface("test_s");
  t2->SetModel(unified);
  t2->SetType(dielectric_dielectric);
  t2->SetFinish(ground);
  t2->SetSigmaAlpha(0.1);
  
  G4OpticalSurface* t3 = new G4OpticalSurface("test_s2");
  t3->SetModel(unified);
  t3->SetType(dielectric_dielectric);
  t3->SetFinish(polished);

  // *
  // Setting the surface type variable for between the 
  //  aluminium reflector and the interim vacuum.
  //  Using an external file in the /src directory because
  //  GEANT4 macros are ridiculous.
  // *
  char label = readExtFile("surface");
  G4OpticalSurface* sides; 
  G4OpticalSurface* tb;
  
  if(label == '3'){ // POLISHED
	sides = t3; tb = t3;

	G4cout << "POLISHED" << G4endl;

  } else if (label == '2'){ // GROUND
	sides = t2; tb = t2;

	G4cout << "GROUND" << G4endl;

  } else { // BOTH
	sides = t3; tb = t2;

	G4cout << "BOTH" << G4endl;
  }

  // *
  //	TODO: property of the interum vacuum sides for below?
  // *

  // *
  // Set the surface properties on both sides of the interface
  //  between the interim vacuum and the LMO crystal
  // *
  new G4LogicalBorderSurface("Surface2",IV1,fLMO1,sides);
  new G4LogicalBorderSurface("Surface3",IV2,fLMO1,sides);
  new G4LogicalBorderSurface("Surface4",IV3,fLMO1,sides);
  new G4LogicalBorderSurface("Surface5",IV4,fLMO1,sides);
  new G4LogicalBorderSurface("Surface6",IV5,fLMO1,tb);
  new G4LogicalBorderSurface("Surface7",IV6,fLMO1,tb);
 
  new G4LogicalBorderSurface("Surface2",fLMO1,IV1,sides);
  new G4LogicalBorderSurface("Surface3",fLMO1,IV2,sides);
  new G4LogicalBorderSurface("Surface4",fLMO1,IV3,sides);
  new G4LogicalBorderSurface("Surface5",fLMO1,IV4,sides);
  new G4LogicalBorderSurface("Surface6",fLMO1,IV5,tb);
  new G4LogicalBorderSurface("Surface7",fLMO1,IV6,tb);

  // *
  // Set the surface properties on both sides of the interface
  //  between the interim vacuum and the aluminium reflector 
  //  (so the aluminium reflector side is always polished)
  // *
  new G4LogicalBorderSurface("Surface2",AR1,IV1,t3);
  new G4LogicalBorderSurface("Surface3",AR2,IV2,t3);
  new G4LogicalBorderSurface("Surface4",AR3,IV3,t3);
  new G4LogicalBorderSurface("Surface5",AR4,IV4,t3);
  new G4LogicalBorderSurface("Surface6",LD5,IV5,t3);
  new G4LogicalBorderSurface("Surface7",LD6,IV6,t3);
 
  new G4LogicalBorderSurface("Surface2",IV1,AR1,t3);
  new G4LogicalBorderSurface("Surface3",IV2,AR2,t3);
  new G4LogicalBorderSurface("Surface4",IV3,AR3,t3);
  new G4LogicalBorderSurface("Surface5",IV4,AR4,t3);
  new G4LogicalBorderSurface("Surface6",IV5,LD5,t3);
  new G4LogicalBorderSurface("Surface7",IV6,LD6,t3);
 
  G4cout << G4endl << "Setting Optical Surface at Junctions." << G4endl;



  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  LightDetectorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
  auto fLDSD
    = new LightDetectorSD("LightDetectorSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);
 
  // *
  // Make the Light Detectors and Aluminium Reflectors Sensitive Detectors 
  // *
  SetSensitiveDetector("LightDetectorLV",fLDSD);
  SetSensitiveDetector("ARLV",fLDSD);

  // *
  // Now we can set the property of being a sensitive detector for the LMO 
  //  with an external file (gosh darn macros are a pain to figure out)
  // *
  char label = readExtFile("sensitive");
  
  if(label == '1'){ 	
	SetSensitiveDetector("LMO",fLDSD);
	SetSensitiveDetector("LMO_primary",fLDSD);

	G4cout << "LMO sensitive" << G4endl;
  } else {
	G4cout << "LMO not sensitive" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && WorldMaterial != pmat)
  {
    WorldMaterial = pmat;
    if(worldLV)
    {
      worldLV->SetMaterial(WorldMaterial);
      WorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "World material set to " << WorldMaterial->GetName() << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetLMOMaterial(const G4String& mat)
{
  G4cout << "Attempting to set LMO material to " << LMOMaterial->GetName() << G4endl;
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && LMOMaterial != pmat)
  {
    LMOMaterial = pmat;
    if(fLMO_LV)
    {
      fLMO_LV->SetMaterial(LMOMaterial);
      LMOMaterial->SetMaterialPropertiesTable(fLMOMPT);
      // LMOMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "LMO material set to " << LMOMaterial->GetName() << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetLDMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && LightDetectorMaterial != pmat)
  {
    LightDetectorMaterial = pmat;
    if(LightDetectorLV)
    {
      LightDetectorLV->SetMaterial(LightDetectorMaterial);
      LightDetectorMaterial->SetMaterialPropertiesTable(fLDMPT);
      // LMOMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "LMO material set to " << LightDetectorMaterial->GetName() << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const G4String& prop,
                                       G4MaterialPropertyVector* mpv)
{
  fWorldMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddLMOMPV(const G4String& prop,
                                      G4MaterialPropertyVector* mpv)
{
  fLMOMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the LMO is now: " << G4endl;
  fLMOMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddLDMPV(const G4String& prop,
                                      G4MaterialPropertyVector* mpv)
{
  fLDMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the LD is now: " << G4endl;
  fLDMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPC(const G4String& prop, G4double v)
{
  fWorldMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddLMOMPC(const G4String& prop, G4double v)
{
  fLMOMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the LMO is now: " << G4endl;
  fLMOMPT->DumpTable();
  G4cout << "............." << G4endl;
}
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// void DetectorConstruction::AddLDMPC(const G4String& prop, G4double v)
// {
//   fLDMPT->AddConstProperty(prop, v);
//   G4cout << "The MPT for the LD is now: " << G4endl;
//   fLDMPT->DumpTable();
//   G4cout << "............." << G4endl;
// }










// *
// Dude, GEANT4 macros are impossible to figure out ........ so use an external file
//  in the /src folder for that logick (essentially the same thing except that macros
//  are more obtuse and more difficult than like 5 LoC >:( )
// *
char readExtFile(const char* filename){

	FILE *p = fopen(filename,"r");
	
	if(p == NULL){
		return '0';
	}

	return fgetc(p);
}






















