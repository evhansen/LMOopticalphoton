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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
// G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

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

  G4double worldabslength[3] = {100*mm, 100*mm, 100*mm};
  G4double energies[3] = {2.0*eV, 3.0*eV, 4.0*eV};

  fLMOMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
  fLMOMPT->AddProperty("RINDEX", energies, rindex, numentries);

  fLDMPT->AddProperty("RINDEX", energies, LDrindex, numentries);
  fLDMPT->AddProperty("ABSLENGTH", energies, LDabslength, numentries);

  // fWorldMPT->AddProperty("ABSLENGTH",energies, worldabslength, numentries);

  fSurface = new G4OpticalSurface("Surface");
  fSurface->SetType(dielectric_dielectric);
  fSurface->SetFinish(polished);
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

  if ( ! WorldMaterial || ! LightDetectorMaterial || ! LMOMaterial ) {
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
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // LMO
  //
  auto LMO_box
    = new G4Box("LMO",          // its name
                fLMO_xy, fLMO_xy, fLMO_xy);    //its size

  fLMO_LV
    = new G4LogicalVolume(
                LMO_box,         // its solid
                LMOMaterial,    // its material
                "LMO");          // its name


  // G4PVPlacement* fLMO1 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO1",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO2 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,-fLMO_C2C_inTower,0),
                fLMO_LV,          // its logical volume
                "LMO2",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                1,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO3 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,0,0),
                fLMO_LV,          // its logical volume
                "LMO3",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                2,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO4 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,fLMO_C2C_diffTower,0),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO4",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                3,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO5 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,0,0),
                fLMO_LV,          // its logical volume
                "LMO5",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                4,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO6 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,0),
                fLMO_LV,          // its logical volume
                "LMO6",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                5,                // copy number
                fCheckOverlaps);  // checking overlaps


  // G4PVPlacement* fLMO7 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,0),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO7",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                6,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO8 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,0),
                fLMO_LV,          // its logical volume
                "LMO8",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                7,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO9 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,0),
                fLMO_LV,          // its logical volume
                "LMO9",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                8,                // copy number
                fCheckOverlaps);  // checking overlaps


  // G4PVPlacement* fLMO10 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,0,fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO10",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                9,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO11 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,-fLMO_C2C_inTower,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO11",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                10,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO12 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,0,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO12",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                11,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO13 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,fLMO_C2C_diffTower,fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO13",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                12,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO14 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,0,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO14",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                13,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO15 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO15",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                14,                // copy number
                fCheckOverlaps);  // checking overlaps


  // G4PVPlacement* fLMO16 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO16",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                15,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO17 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO17",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                16,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO18 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO18",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                17,                // copy number
                fCheckOverlaps);  // checking overlaps



  // G4PVPlacement* fLMO10 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,0,-fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO19",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                18,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO11 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,-fLMO_C2C_inTower,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO20",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                19,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO12 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,0,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO21",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                20,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO13 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,fLMO_C2C_diffTower,-fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO22",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                21,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO14 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,0,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO23",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                22,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO15 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO24",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                23,                // copy number
                fCheckOverlaps);  // checking overlaps


  // G4PVPlacement* fLMO16 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,-fLMO_UpperRowZ),  // place off center eventually
                fLMO_LV,          // its logical volume
                "LMO25",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                24,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO17 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO26",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                25,                // copy number
                fCheckOverlaps);  // checking overlaps

  // G4PVPlacement* fLMO18 =
    new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,-fLMO_UpperRowZ),
                fLMO_LV,          // its logical volume
                "LMO27",          // its name
                worldLV,                // its mother  volume
                false,            // no boolean operation
                26,                // copy number
                fCheckOverlaps);  // checking overlaps
  //
  // LightDetector
  //
  auto LightDetector
    = new G4Box("LightDetector",     // its name
                 fLD_xy, fLD_xy, fLD_z); // its size

  LightDetectorLV
    = new G4LogicalVolume(
                 LightDetector,     // its solid
                 LightDetectorMaterial,  // its material
                 "LightDetectorLV");   // its name

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,fLD_C2C_Top),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD1T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               0,                // copy number
               fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,-fLD_C2C_Bot),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD1B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               1,                // copy number
               fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,-fLMO_C2C_inTower,fLD_C2C_Top),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD2T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              2,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,-fLMO_C2C_inTower,-fLD_C2C_Bot),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD2B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              3,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,0,fLD_C2C_Top),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD3T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               4,                // copy number
               fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,0,-fLD_C2C_Bot),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD3B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               5,                // copy number
               fCheckOverlaps);  // checking overlaps

 new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,fLMO_C2C_diffTower,fLD_C2C_Top),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD4T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              6,                // copy number
              fCheckOverlaps);  // checking overlaps

 new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,fLMO_C2C_diffTower,-fLD_C2C_Bot),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD4B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              7,                // copy number
              fCheckOverlaps);  // checking overlaps


 new G4PVPlacement(
             0,                // no rotation
             G4ThreeVector(-fLMO_C2C_diffTower,0,fLD_C2C_Top),  // at (0,0,0)
             LightDetectorLV,          // its logical volume
             "LD5T",    // its name
             worldLV,          // its mother  volume
             false,            // no boolean operation
             8,                // copy number
             fCheckOverlaps);  // checking overlaps

 new G4PVPlacement(
             0,                // no rotation
             G4ThreeVector(-fLMO_C2C_diffTower,0,-fLD_C2C_Bot),  // at (0,0,0)
             LightDetectorLV,          // its logical volume
             "LD5B",    // its name
             worldLV,          // its mother  volume
             false,            // no boolean operation
             9,                // copy number
             fCheckOverlaps);  // checking overlaps

 new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,fLD_C2C_Top),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD6T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              10,                // copy number
              fCheckOverlaps);  // checking overlaps

 new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,-fLD_C2C_Bot),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD6B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              11,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,fLD_C2C_Top),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD7T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               12,                // copy number
               fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,-fLD_C2C_Bot),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD7B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               13,                // copy number
               fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,fLD_C2C_Top),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD8T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              14,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,-fLD_C2C_Bot),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD8B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              15,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,fLD_C2C_Top),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD9T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               16,                // copy number
               fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,-fLD_C2C_Bot),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD9B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               17,                // copy number
               fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,0,fLD_UpperRowZ),  // at (0,0,0)
                LightDetectorLV,          // its logical volume
                "LD10T",    // its name
                worldLV,          // its mother  volume
                false,            // no boolean operation
                18,                // copy number
                fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,-fLMO_C2C_inTower,fLD_UpperRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD11T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               19,                // copy number
               fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,0,fLD_UpperRowZ),  // at (0,0,0)
                LightDetectorLV,          // its logical volume
                "LD12T",    // its name
                worldLV,          // its mother  volume
                false,            // no boolean operation
                20,                // copy number
                fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,fLMO_C2C_diffTower,fLD_UpperRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD13T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               21,                // copy number
               fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(-fLMO_C2C_diffTower,0,fLD_UpperRowZ),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD14T",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              22,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,fLD_UpperRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD15T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               23,                // copy number
               fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,fLD_UpperRowZ),  // at (0,0,0)
                LightDetectorLV,          // its logical volume
                "LD16T",    // its name
                worldLV,          // its mother  volume
                false,            // no boolean operation
                24,                // copy number
                fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,fLD_UpperRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD17T",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               25,                // copy number
               fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,fLD_UpperRowZ),  // at (0,0,0)
                LightDetectorLV,          // its logical volume
                "LD18T",    // its name
                worldLV,          // its mother  volume
                false,            // no boolean operation
                26,                // copy number
                fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,-fLD_LowerRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD1B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               1,                // copy number
               fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,-fLMO_C2C_inTower,-fLD_LowerRowZ),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD2B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              3,                // copy number
              fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,0,-fLD_LowerRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD3B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               5,                // copy number
               fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(0,fLMO_C2C_diffTower,-fLD_LowerRowZ),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD4B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              7,                // copy number
              fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
             0,                // no rotation
             G4ThreeVector(-fLMO_C2C_diffTower,0,-fLD_LowerRowZ),  // at (0,0,0)
             LightDetectorLV,          // its logical volume
             "LD5B",    // its name
             worldLV,          // its mother  volume
             false,            // no boolean operation
             9,                // copy number
             fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(-fLMO_C2C_diffTower,-fLMO_C2C_inTower,-fLD_LowerRowZ),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD6B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              11,                // copy number
              fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(fLMO_C2C_inTower,-fLMO_C2C_inTower,-fLD_LowerRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD7B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               13,                // copy number
               fCheckOverlaps);  // checking overlaps


  new G4PVPlacement(
              0,                // no rotation
              G4ThreeVector(fLMO_C2C_inTower,fLMO_C2C_diffTower,-fLD_LowerRowZ),  // at (0,0,0)
              LightDetectorLV,          // its logical volume
              "LD8B",    // its name
              worldLV,          // its mother  volume
              false,            // no boolean operation
              15,                // copy number
              fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(-fLMO_C2C_diffTower,fLMO_C2C_diffTower,-fLD_LowerRowZ),  // at (0,0,0)
               LightDetectorLV,          // its logical volume
               "LD9B",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               17,                // copy number
               fCheckOverlaps);  // checking overlaps



  // G4LogicalBorderSurface* surface1
  //   = new G4LogicalBorderSurface("Surface", fLMO1, worldPV, fSurface);
  // G4OpticalSurface* opticalSurface1
  //   = dynamic_cast<G4OpticalSurface*>(surface1->GetSurface(fLMO1, worldPV)->GetSurfaceProperty());


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
  SetSensitiveDetector("LightDetectorLV",fLDSD);

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
