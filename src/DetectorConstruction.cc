#include "DetectorConstruction.hh"
#include "LightDetectorSD.hh"
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
	
	DefineMaterials();
	
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
	G4VPhysicalVolume* wpv = DefineWorld();
	DefineVolumes();
	return wpv;
}



void DetectorConstruction::DefineMaterials()
{
	char ldmat = readExtFile("./options/ldmat");
	char emat = readExtFile("./options/emat");
	char crysmat = readExtFile("./options/crysmat");
	//char arcmat = readExtFile("./options/arcmat");
	
	// optical surface MPT (because of k -.-)
	LDk = new G4MaterialPropertiesTable();
	EMk = new G4MaterialPropertiesTable();
	ARCk = new G4MaterialPropertiesTable();


	//************
	// Define Materials 
	//************
	
	DefineWorldM();	

	DefineCrysM(crysmat);

	DefineLDM(ldmat);

	DefineEMM(emat);

	//DefineARCM();


	// *****************
	// TODO: Optical surface MPTs?
	// *****************


	G4cout << *(G4Material::GetMaterialTable()) << G4endl;


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
				900,                // copy number
				fCheckOverlaps);  // checking overlaps


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	return worldPV;
}



void DetectorConstruction::DefineVolumes()
{	
	float ldg = stof(readExtFileLine("./options/ldg"));
	float emg = stof(readExtFileLine("./options/emg")); 

	float arct = stof(readExtFileLine("./options/arct"));
	//char arc = readExtFile("./options/arc");
	

	// *********
	// Defining Sizes and Pos' 
	// *********


	// TODO: throw in crysd! (pass to PAG)
	//G4double LMOHalfSizes[3]
	//	= {1*cm, 1*cm, 1*cm}; // 10 mm x 10 mm x 10 mm
	G4double LMOHalfSizes[3]
		= {22.5*mm, 22.5*mm, 22.5*mm};

	fLMO_xy = LMOHalfSizes[0]; // what if [0] != [1]? TODO: get rid of this

	// note: cannot change LMOHS' then expect below to work + don't change units
	G4double LDHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double EMHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double ARCHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], arct*mm};


	// not necessary but it's ok
	G4double LDgap, EMgap;
	LDgap = ((G4double)ldg)*mm;
	EMgap = ((G4double)emg)*mm;


	G4double LDPos[19]
		= {0, 0, LMOHalfSizes[2] + LDgap + LDHalfSizes[2],
			0, 0, -LMOHalfSizes[2] - LDgap - LDHalfSizes[2],
			0, -LMOHalfSizes[1] - LDgap - LDHalfSizes[2], 0,
			0, LMOHalfSizes[1] + LDgap + LDHalfSizes[2], 0,
			-LMOHalfSizes[0] - LDgap - LDHalfSizes[2], 0, 0,
			LMOHalfSizes[0] + LDgap + LDHalfSizes[2], 0, 0}; 
	// All LDHS' because rotated appropriately


	G4double EMPos[19]
		= {0, 0, LMOHalfSizes[2] + EMgap + EMHalfSizes[2],
			0, 0, -LMOHalfSizes[2] - EMgap - EMHalfSizes[2],
			0, -LMOHalfSizes[1] - EMgap - EMHalfSizes[2], 0,
			0, LMOHalfSizes[1] + EMgap + EMHalfSizes[2], 0,
			-LMOHalfSizes[0] - EMgap - EMHalfSizes[2], 0, 0,
			LMOHalfSizes[0] + EMgap + EMHalfSizes[2], 0, 0};


	// TODO: polish ARC after particles are fixed
	/*if(arc == '1'){
		LDPos[2] += 2*ARCHalfSizes[2];
		LDPos[5] -= 2*ARCHalfSizes[2];
		LDPos[7] -= 2*ARCHalfSizes[2];
		LDPos[10] += 2*ARCHalfSizes[2];
		LDPos[12] -= 2*ARCHalfSizes[2];
		LDPos[15] += 2*ARCHalfSizes[2];

		EMPos[2] += 2*ARCHalfSizes[2];
		EMPos[5] -= 2*ARCHalfSizes[2];
		EMPos[7] -= 2*ARCHalfSizes[2];
		EMPos[10] += 2*ARCHalfSizes[2];
		EMPos[12] -= 2*ARCHalfSizes[2];
		EMPos[15] += 2*ARCHalfSizes[2];
	}*/

	// TODO: determine which slots need EMgap 
	/*G4double ARCPos[19]
		= {0, 0, LMOHalfSizes[2] + LDgap + ARCHalfSizes[2],
			0, 0, -LMOHalfSizes[2] - LDgap - ARCHalfSizes[2],
			0, -LMOHalfSizes[1] - LDgap - ARCHalfSizes[2], 0,
			0, LMOHalfSizes[1] + LDgap + ARCHalfSizes[2], 0,
			-LMOHalfSizes[0] - LDgap - ARCHalfSizes[2], 0, 0,
			LMOHalfSizes[0] + LDgap + ARCHalfSizes[2], 0, 0};
	*/



	G4double IVHalfSizes[4]; // [2],[3] for LDgap and EMgap respectively
	G4double IVPos[19] = {0};

	IVHalfSizes[0] = LMOHalfSizes[0]; // x 
	IVHalfSizes[1] = LMOHalfSizes[1]; // y

	IVHalfSizes[2] = (LDgap / 2); // z
	IVHalfSizes[3] = (EMgap / 2); // z

	
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");

	
	// Yeah, yeah, refactored later
	if(face1 == 'l'){ //face1 has an LD
		IVPos[2] = (LDgap / 2) + LMOHalfSizes[2]; 

	} else if(face1 == 'e'){ //face1 has an EM	
		IVPos[2] = (EMgap / 2) + LMOHalfSizes[2]; 

	} else { //face 1 has nothing
		IVPos[2] = LMOHalfSizes[2]; 
	}

	if(face2 == 'l'){ //LD
		IVPos[5] = -(LDgap / 2) - LMOHalfSizes[2];

	} else if(face2 == 'e'){ // EM
		IVPos[5] = -(EMgap / 2) - LMOHalfSizes[2];

	} else { //nothing
		IVPos[5] = - LMOHalfSizes[2];
	}

	if(face3 == 'l'){ //LD
		IVPos[7] = -(LDgap / 2) - LMOHalfSizes[1]; 

	} else if(face3 == 'e'){ // EM
		IVPos[7] = -(EMgap / 2) - LMOHalfSizes[1]; 

	} else { //nothing
		IVPos[7] = - LMOHalfSizes[1]; 
	}

	if(face4 == 'l'){ //LD
		IVPos[10] = (LDgap / 2) + LMOHalfSizes[1];

	} else if(face4 == 'e'){ // EM
		IVPos[10] = (EMgap / 2) + LMOHalfSizes[1]; 

	} else {//nothing
		IVPos[10] = LMOHalfSizes[1]; 
	}

	if(face5 == 'l'){ //LD
		IVPos[12] = -(LDgap / 2) - LMOHalfSizes[0]; 

	} else if(face4 == 'e'){ // EM
		IVPos[12] = -(EMgap / 2) - LMOHalfSizes[0];

	} else { //nothing
		IVPos[12] = - LMOHalfSizes[0];
	}

	if(face6 == 'l'){ //LD
		IVPos[15] = (LDgap / 2) + LMOHalfSizes[0];

	} else if(face4 == 'e'){ // EM
		IVPos[15] = (EMgap / 2) + LMOHalfSizes[0];

	} else { //nothing
		IVPos[15] = LMOHalfSizes[0];
	}



	// *****************
	// Defining Volumes
	// *****************


	G4VPhysicalVolume* fLMO1 = DefineCrys(LMOHalfSizes);
	

	G4VPhysicalVolume* IV1,*IV2,*IV3,*IV4,*IV5,*IV6;
	
	DefineIV(IVHalfSizes,IVPos,
				&IV1,&IV2,&IV3,&IV4,&IV5,&IV6);


	G4VPhysicalVolume* LD1,*LD2,*LD3,*LD4,*LD5,*LD6;
	
	DefineLD(LDHalfSizes,LDPos,
				&LD1,&LD2,&LD3,&LD4,&LD5,&LD6);

		
	G4VPhysicalVolume* EM1,*EM2,*EM3,*EM4,*EM5,*EM6;
	
	DefineEM(EMHalfSizes,EMPos,
				&EM1,&EM2,&EM3,&EM4,&EM5,&EM6);


	//G4VPhysicalVolume* ARC1,*ARC2,*ARC3,*ARC4,*ARC5,*ARC6;

	//DefineARC(ARCHalfSizes,ARCPos,
	//		&ARC1,&ARC2,&ARC3,&ARC4,&ARC5,&ARC6);	

	
	// *****************
	// Defining Surfaces
	// *****************

	//TODO: add back surfaces when particles are fixed

}




void DetectorConstruction::ConstructSDandField()
{
	//TODO: I know I know cringe
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");
	char lmos = readExtFile("./options/lmos");

	
	G4VSensitiveDetector* fLDSD 
		= new LightDetectorSD("LightDetectorSD", 
				"AbsorberHitsCollection", fNofLayers);

	G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);


	if(lmos == 's'){
		SetSensitiveDetector("LMO_primary",fLDSD);
	
		G4cout << "LMO sensitive" << G4endl;
	}
	
	
	if(face1 == 'l' || face2 == 'l' || face3 == 'l' || face4 == 'l' || face5 == 'l' || face6 == 'l'){
		SetSensitiveDetector("LDLV",fLDSD); // LDs sensitive
	
		G4cout << "LDs sensitive" << G4endl;
	}

	
	//not really necessary
	//SetSensitiveDetector("EMLV",fLDSD);
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

