#include "DetectorConstruction.hh"
#include "LightDetectorSD.hh"
#include "ELightDetectorSD.hh"
#include "DetectorMessenger.hh"

#include "G4GeometryManager.hh"
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
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4UserLimits.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "DefMats.hh"
#include "DefGeos.hh"


#define a 1

G4ThreadLocal


static char readExtFile(const char* filename);
std::string readExtFileLine(const char* filename);


DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fCheckOverlaps(true),
fDetectorMessenger(nullptr),
fNofLayers(-1)
{
	fNofLayers = 1;
	
	
	fDetectorMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{
	delete worldLV;
	
	delete WorldMaterial;
	delete LMOMaterial;
	delete LightDetectorMaterial;
	delete EMMaterial;
	delete ARCMaterial;

	// MPTs for surfaces
	delete LDk;
	delete EMk;
	delete ARCk;
	
	delete fDetectorMessenger;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

	DefineMaterials();
	G4VPhysicalVolume* wpv = DefineWorld();
	DefineVolumes();
	return wpv;
}



void DetectorConstruction::DefineMaterials()
{
	char ldmat = readExtFile("./options/ldmat");
	char emat = readExtFile("./options/emat");
	char crysmat = readExtFile("./options/crysmat");
	char arc = readExtFile("./options/arcmat");

#if a
	G4cout << "\nLD material: " << ldmat 
	<< "\nEM material: " << emat 
	<< "\nCrys material: " << crysmat << "\n" << G4endl;
#endif
	
	// optical surface MPT (because of k -.-)
	LDk = new G4MaterialPropertiesTable();
	EMk = new G4MaterialPropertiesTable();
	ARCk = new G4MaterialPropertiesTable();


	//************
	// Define Materials 
	//************

#if a
	G4cout << "\nDefining world materials.\n" << G4endl;
#endif
	DefineWorldM();	

#if a
	G4cout << "\nDefining crys materials.\n" << G4endl;
#endif
	DefineCrysM(crysmat);

#if a
	G4cout << "\nDefining LD materials.\n" << G4endl;
#endif
	DefineLDM(ldmat);

#if a
	G4cout << "\nDefining EM materials.\n" << G4endl;
#endif
	DefineEMM(emat);


	if(arc == '1' || arc == '2' || arc == 'n'){
#if a
		G4cout << "\nDefining ARC materials.\n" << G4endl;
#endif
		DefineARCM(arc);
	}


	// *****************
	// TODO: Optical surface MPTs?
	// *****************


	// segfaults on *wo4?
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}




G4VPhysicalVolume* DetectorConstruction::DefineWorld(){

	G4double worldSizeXY = 10. * m;

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldSizeXY);

	auto worldS
		= new G4Box("World",           // its name
				worldSizeXY/2, worldSizeXY/2, worldSizeXY/2); // its size

	worldLV
		= new G4LogicalVolume(
				worldS,           // its solid
				WorldMaterial,  // its material
				"World");         // its name

	G4VPhysicalVolume* worldPV
		= new G4PVPlacement(
				0,                // no rotation
				G4ThreeVector(),  // at (0,0,0)
				worldLV,          // its logical volume
				"World",          // its name
				0,                // its mother  volume
				false,            // no boolean operation
				100,                // copy number
				fCheckOverlaps);  // checking overlaps


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	return worldPV;
}



void DetectorConstruction::DefineVolumes()
{	
	float ldg = stof(readExtFileLine("./options/ldg"));
	float emg = stof(readExtFileLine("./options/emg")); 

	float arct = stof(readExtFileLine("./options/arct"));
	char arc = readExtFile("./options/arcmat");


	G4Material* CuMaterial;
	G4Material* AlMaterial;
	DefineMisc(&CuMaterial, &AlMaterial);

	G4double LMOHalfSizes[3] = {1*cm, 1*cm, 1*cm};
	fLMO_xy = LMOHalfSizes[0]; // what if [0] != [1]? TODO: get rid of this
	
	G4VPhysicalVolume* fLMO1 = DefineCrys(LMOHalfSizes);

	G4VisAttributes* IVVA = new G4VisAttributes();
	IVVA->SetColor(1.0,0.0,1.0,0.15);
	IVVA->SetForceSolid(true);

	G4double LDzPos = LMOHalfSizes[2] + 0.2*cm + (0.028/2)*cm;

	G4int wR = 2.54*cm;
	
	// CYLINDRICAL
	G4VSolid* LDwafer = new G4Tubs("LDwafer",0.*cm,wR,0.14*mm,0.*rad,2.*pi*rad);
	G4LogicalVolume* LDLV = new G4LogicalVolume(LDwafer,LightDetectorMaterial,"LDLV");
	LDLV->SetVisAttributes(IVVA);

	G4VisAttributes* IVVA2 = new G4VisAttributes();
	IVVA2->SetColor(1,0,0,0.15);
	IVVA2->SetForceSolid(true);

	G4LogicalVolume* ELDLV = new G4LogicalVolume(LDwafer,LightDetectorMaterial,"ELDLV");
	ELDLV->SetVisAttributes(IVVA2);


	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);

	//new G4PVPlacement(t1,G4ThreeVector(0,LDzPos,0),LDLV,"LD1",worldLV,false,10,fCheckOverlaps);	

	new G4PVPlacement(t1,G4ThreeVector(0,LDzPos,0),LDLV,"LD3",worldLV,false,13,fCheckOverlaps);	


	new G4PVPlacement(0,G4ThreeVector(0,0,LDzPos),ELDLV,"LD1",worldLV,false,10,fCheckOverlaps);	
	new G4PVPlacement(t1,G4ThreeVector(0,-LDzPos,0),ELDLV,"LD2",worldLV,false,12,fCheckOverlaps);	
	new G4PVPlacement(p1,G4ThreeVector(-LDzPos,0,0),ELDLV,"LD4",worldLV,false,14,fCheckOverlaps);	
	new G4PVPlacement(p1,G4ThreeVector(LDzPos,0,0),ELDLV,"LD5",worldLV,false,15,fCheckOverlaps);	


	G4VSolid* XInterimVacuum = new G4Box("LDInterimVacuum",LMOHalfSizes[0],LMOHalfSizes[1],(0.2/2)*cm); 

	G4LogicalVolume* LDIVLV = new G4LogicalVolume(XInterimVacuum,WorldMaterial,"XInterimVacuumLV"); 

	G4VisAttributes* IVV = new G4VisAttributes();
	IVV->SetColor(1.0,0.0,1.0,0.15);
	IVV->SetForceSolid(true);
	LDIVLV->SetVisAttributes(IVV);

	G4double IVPos[19]
		= {0, 0,LMOHalfSizes[0] + 0.1*cm,
		0, 0, -LMOHalfSizes[0] - 0.1*cm,
		0, -LMOHalfSizes[0] - 0.1*cm, 0,
		0, LMOHalfSizes[0] + 0.1*cm, 0,
		-LMOHalfSizes[0] - 0.1*cm, 0, 0,
		LMOHalfSizes[0] + 0.1*cm, 0, 0}; 

/*	new G4PVPlacement(0,G4ThreeVector(IVPos[0],IVPos[1],IVPos[2]),LDIVLV,"IV1T",worldLV,false,50,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(IVPos[3],IVPos[4],IVPos[5]),LDIVLV,"IV1T",worldLV,false,51,fCheckOverlaps);
	new G4PVPlacement(t1,G4ThreeVector(IVPos[6],IVPos[7],IVPos[8]),LDIVLV,"IV1T",worldLV,false,52,fCheckOverlaps);
	new G4PVPlacement(t1,G4ThreeVector(IVPos[9],IVPos[10],IVPos[11]),LDIVLV,"IV1T",worldLV,false,53,fCheckOverlaps);
	new G4PVPlacement(p1,G4ThreeVector(IVPos[12],IVPos[13],IVPos[14]),LDIVLV,"IV1T",worldLV,false,54,fCheckOverlaps);
	new G4PVPlacement(p1,G4ThreeVector(IVPos[15],IVPos[16],IVPos[17]),LDIVLV,"IV1T",worldLV,false,55,fCheckOverlaps);
*/

	G4int dc = 15*cm + 0.05*cm;

	G4double B1Pos[19]
		= {0, 0, wR + dc,
		0, 0, -wR - dc,
		0, -wR - dc, 0,
		0, wR + dc, 0,
		-wR - dc, 0, 0,
		wR + dc, 0, 0}; 

	G4int dc2 = 0.5*cm + (0.05*3)*cm;

	G4double B2Pos[19]
		= {0, 0, B1Pos[2] + dc2,
		0, 0, B1Pos[5] - dc2,
		0, B1Pos[7] - dc2, 0,
		0, B1Pos[10] + dc2, 0,
		B1Pos[12] - dc2, 0, 0,
		B1Pos[15] + dc2, 0, 0}; 

	auto AlBox = new G4Box("AlBox",(2.54+15)*cm,(2.54+15)*cm,0.5*mm);    
	auto CuBox = new G4Box("CuBox",(2.54+15+0.5+(0.05*2))*cm,(2.54+15+0.5+(0.05*2))*cm,0.5*mm);    

	G4VisAttributes* ELDVA = new G4VisAttributes();
	ELDVA->SetColor(0.0,0.5,0.5,0.15);
	ELDVA->SetForceSolid(true);

	G4VisAttributes* ELDVA2 = new G4VisAttributes();
	ELDVA2->SetColor(0.5,0.5,0.0,0.15);
	ELDVA2->SetForceSolid(true);

	G4LogicalVolume* AlBLV = new G4LogicalVolume(AlBox,AlMaterial,"AlBLV"); 
	G4LogicalVolume* CuBLV = new G4LogicalVolume(CuBox,CuMaterial,"CuBLV"); 

	AlBLV->SetVisAttributes(ELDVA);
	CuBLV->SetVisAttributes(ELDVA2);

	new G4PVPlacement(0,G4ThreeVector(B1Pos[0],B1Pos[1],B1Pos[2]),AlBLV,"AlB1",worldLV,false,30,fCheckOverlaps);	

	new G4PVPlacement(0,G4ThreeVector(B1Pos[3],B1Pos[4],B1Pos[5]),AlBLV,"AlB2",worldLV,false,31,fCheckOverlaps);	

	new G4PVPlacement(t1,G4ThreeVector(B1Pos[6],B1Pos[7],B1Pos[8]),AlBLV,"AlB3",worldLV,false,32,fCheckOverlaps);	

	new G4PVPlacement(t1,G4ThreeVector(B1Pos[9],B1Pos[10],B1Pos[11]),AlBLV,"AlB4",worldLV,false,33,fCheckOverlaps);	

	new G4PVPlacement(p1,G4ThreeVector(B1Pos[12],B1Pos[13],B1Pos[14]),AlBLV,"AlB5",worldLV,false,34,fCheckOverlaps);	
	
	new G4PVPlacement(p1,G4ThreeVector(B1Pos[15],B1Pos[16],B1Pos[17]),AlBLV,"AlB6",worldLV,false,35,fCheckOverlaps);	


	new G4PVPlacement(0,G4ThreeVector(B2Pos[0],B2Pos[1],B2Pos[2]),CuBLV,"CuB1",worldLV,false,40,fCheckOverlaps);	

	new G4PVPlacement(0,G4ThreeVector(B2Pos[3],B2Pos[4],B2Pos[5]),CuBLV,"CuB2",worldLV,false,41,fCheckOverlaps);	

	new G4PVPlacement(t1,G4ThreeVector(B2Pos[6],B2Pos[7],B2Pos[8]),CuBLV,"CuB3",worldLV,false,42,fCheckOverlaps);	

	new G4PVPlacement(t1,G4ThreeVector(B2Pos[9],B2Pos[10],B2Pos[11]),CuBLV,"CuB4",worldLV,false,43,fCheckOverlaps);	

	new G4PVPlacement(p1,G4ThreeVector(B2Pos[12],B2Pos[13],B2Pos[14]),CuBLV,"CuB5",worldLV,false,44,fCheckOverlaps);	
	
	new G4PVPlacement(p1,G4ThreeVector(B2Pos[15],B2Pos[16],B2Pos[17]),CuBLV,"CuB6",worldLV,false,45,fCheckOverlaps);	



}




void DetectorConstruction::ConstructSDandField()
{
	G4VSensitiveDetector* fLDSD 
		= new LightDetectorSD("LightDetectorSD", 
				"AbsorberHitsCollection", fNofLayers);

	G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);

	G4VSensitiveDetector* fELDSD 
		= new ELightDetectorSD("ELightDetectorSD", 
				"HitsCollection", fNofLayers);
	
	G4SDManager::GetSDMpointer()->AddNewDetector(fELDSD);
	
	SetSensitiveDetector("LMO_primary",fLDSD);

	SetSensitiveDetector("LDLV",fLDSD);
	SetSensitiveDetector("ELDLV",fELDSD);

	SetSensitiveDetector("AlBLV",fLDSD);
	SetSensitiveDetector("CuBLV",fLDSD);
}





char readExtFile(const char* filename){
	FILE *p = fopen(filename,"r");
	
	if(p == NULL){
		return '0';
	}
	
	char out = fgetc(p);

	fclose(p);

	return out;

}

std::string readExtFileLine(const char* filename){
	std::string out;
	std::ifstream p;

	p.open(filename);

	std::getline(p,out);

	return out;
}

