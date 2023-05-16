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

#define RegularGeometry 0
#define ARC 0
	#define LogicalBorders 0
#define ELD 0

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
	

	// *********
	// Defining Sizes and Pos' 
	// *********



#if RegularGeometry
	// TODO: throw in crysd! (pass to PAG)
	//G4double LMOHalfSizes[3]
	//	= {1*cm, 1*cm, 1*cm}; // 10 mm x 10 mm x 10 mm
	G4double LMOHalfSizes[3]
		//= {10*mm, 10*mm, 10*mm};
		//= {22.5*mm, 22.5*mm, 22.5*mm};
		= {10*cm, 10*cm, 10*cm};

	fLMO_xy = LMOHalfSizes[0]; // what if [0] != [1]? TODO: get rid of this

	
	// note: cannot change LMOHS' then expect below to work + don't change units
	G4double LDHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double EMHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double ARCHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], arct*nm};


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
	if(arc == '1' || arc == '2' || arc == 'n'){
		LDPos[2] += 2*ARCHalfSizes[2];
		LDPos[5] -= 2*ARCHalfSizes[2];
		LDPos[7] -= 2*ARCHalfSizes[2];
		LDPos[10] += 2*ARCHalfSizes[2];
		LDPos[12] -= 2*ARCHalfSizes[2];
		LDPos[15] += 2*ARCHalfSizes[2];
	}

	// TODO: determine which slots need EMgap 
	G4double ARCPos[19]
		= {0, 0, LMOHalfSizes[2] + LDgap + ARCHalfSizes[2],
			0, 0, -LMOHalfSizes[2] - LDgap - ARCHalfSizes[2],
			0, -LMOHalfSizes[1] - LDgap - ARCHalfSizes[2], 0,
			0, LMOHalfSizes[1] + LDgap + ARCHalfSizes[2], 0,
			-LMOHalfSizes[0] - LDgap - ARCHalfSizes[2], 0, 0,
			LMOHalfSizes[0] + LDgap + ARCHalfSizes[2], 0, 0};
	



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

	if(face2 == 'l'){ 
		IVPos[5] = -(LDgap / 2) - LMOHalfSizes[2];

	} else if(face2 == 'e'){ 
		IVPos[5] = -(EMgap / 2) - LMOHalfSizes[2];

	} else { 
		IVPos[5] = - LMOHalfSizes[2];
	}

	if(face3 == 'l'){ 
		IVPos[7] = -(LDgap / 2) - LMOHalfSizes[1]; 

	} else if(face3 == 'e'){
		IVPos[7] = -(EMgap / 2) - LMOHalfSizes[1]; 

	} else { 
		IVPos[7] = - LMOHalfSizes[1]; 
	}

	if(face4 == 'l'){ 
		IVPos[10] = (LDgap / 2) + LMOHalfSizes[1];

	} else if(face4 == 'e'){ 
		IVPos[10] = (EMgap / 2) + LMOHalfSizes[1]; 

	} else {
		IVPos[10] = LMOHalfSizes[1]; 
	}

	if(face5 == 'l'){ 
		IVPos[12] = -(LDgap / 2) - LMOHalfSizes[0]; 

	} else if(face4 == 'e'){ 
		IVPos[12] = -(EMgap / 2) - LMOHalfSizes[0];

	} else { 
		IVPos[12] = - LMOHalfSizes[0];
	}

	if(face6 == 'l'){ 
		IVPos[15] = (LDgap / 2) + LMOHalfSizes[0];

	} else if(face4 == 'e'){ 
		IVPos[15] = (EMgap / 2) + LMOHalfSizes[0];

	} else { 
		IVPos[15] = LMOHalfSizes[0];
	}



	// *****************
	// Defining Volumes
	// *****************


	#if a
	G4cout << "\nDefining crys volumes.\n" << G4endl;
	#endif
	G4VPhysicalVolume* fLMO1 = DefineCrys(LMOHalfSizes);
	

	G4VPhysicalVolume* IV1,*IV2,*IV3,*IV4,*IV5,*IV6;
	
	#if a
	G4cout << "\nDefining IV volumes.\n" << G4endl;
	#endif
	DefineIV(IVHalfSizes,IVPos,
				&IV1,&IV2,&IV3,&IV4,&IV5,&IV6);


	G4VPhysicalVolume* LD1,*LD2,*LD3,*LD4,*LD5,*LD6;
	
	#if a
	G4cout << "\nDefining LD volumes.\n" << G4endl;
	#endif
	DefineLD(LDHalfSizes,LDPos,
				&LD1,&LD2,&LD3,&LD4,&LD5,&LD6);

		
	G4VPhysicalVolume* EM1,*EM2,*EM3,*EM4,*EM5,*EM6;
	
	#if a
	G4cout << "\nDefining EM volumes.\n" << G4endl;
	#endif
	DefineEM(EMHalfSizes,EMPos,
				&EM1,&EM2,&EM3,&EM4,&EM5,&EM6);

#endif

	// *****************
	// Defining Surfaces
	// *****************

#if ARC 
	if(arc == '1' || arc == '2' || arc == 'n'){

	#if a
		G4cout << "\nDefining ARC volumes.\n" << G4endl;
	#endif
		G4VPhysicalVolume* ARC1,*ARC2,*ARC3,*ARC4,*ARC5,*ARC6;

		DefineARC(ARCHalfSizes,ARCPos,
				&ARC1,&ARC2,&ARC3,&ARC4,&ARC5,&ARC6);	

	#if LogicalBorders 
		double en[5] = {0.005*eV,1*eV,2*eV,3*eV};
		//double trans2[5] = {0.2,0.6,0.5,0.6,0.5};
		//double trans2[5] = {1,1,1,1,1};
		double trans[5] = {0,0,0,0,0};
		double trans2[5] = {0.5,0.5,0.5,0.5,0.5};

		ARCk->AddProperty("TRANSMITTANCE",en,trans,5);
		//ARCk->AddProperty("REFLECTIVITY",en,ref,5);

		G4MaterialPropertiesTable* ARCd = new G4MaterialPropertiesTable();
		//ARCd->AddProperty("TRANSMITTANCE",en,trans2,5);

		G4OpticalSurface* SiO = new G4OpticalSurface("SiO_s");
		SiO->SetModel(unified);
		SiO->SetType(dielectric_dielectric);
		SiO->SetFinish(polished);
		SiO->SetMaterialPropertiesTable(ARCk);

		G4OpticalSurface* SiO2 = new G4OpticalSurface("SiO2_s");
		SiO2->SetModel(unified);
		SiO2->SetType(dielectric_dielectric);
		//SiO2->SetFinish(ground);
		//SiO2->SetSigmaAlpha(0);
		SiO2->SetFinish(polished);
		SiO2->SetMaterialPropertiesTable(ARCd);


		G4OpticalSurface* Gs = new G4OpticalSurface("G_s");
		Gs->SetModel(glisur);
		Gs->SetType(dielectric_dielectric);
		Gs->SetFinish(ground);
		Gs->SetPolish(0.8);


		if(face1 == 'l') { 
			G4cout << "\nFace 1: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC1",IV1,ARC1,SiO);
				new G4LogicalBorderSurface("ARCLD1",ARC1,LD1,SiO);
				//new G4LogicalBorderSurface("ARCIV1",ARC1,IV1,SiO2);
				new G4LogicalBorderSurface("LDARC1",LD1,ARC1,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC1",fLMO1,ARC1,SiO);
				new G4LogicalBorderSurface("ARCLD1",ARC1,LD1,SiO);
				//new G4LogicalBorderSurface("ARCLMO1",ARC1,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC1",LD1,ARC1,SiO2);
			}
		} else if(face1 == 'e'){
			G4cout << "\nFace 1: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM1",fLMO1,EM1,Gs);
		}
	
		if(face2 == 'l'){ 
			G4cout << "\nFace 2: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC2",IV2,ARC2,SiO);
				new G4LogicalBorderSurface("ARCLD2",ARC2,LD2,SiO);
				//new G4LogicalBorderSurface("ARCIV2",ARC2,IV2,SiO2);
				new G4LogicalBorderSurface("LDARC2",LD2,ARC2,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC2",fLMO1,ARC2,SiO);
				new G4LogicalBorderSurface("ARCLD2",ARC2,LD2,SiO);
				//new G4LogicalBorderSurface("ARCLMO2",ARC2,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC2",LD2,ARC2,SiO2);
			}
		} else if(face2 == 'e'){
			G4cout << "\nFace 2: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM2",fLMO1,EM2,Gs);
		}

		if(face3 == 'l'){ 
			G4cout << "\nFace 3: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC3",IV3,ARC3,SiO);
				new G4LogicalBorderSurface("ARCLD3",ARC3,LD3,SiO);
				//new G4LogicalBorderSurface("ARCIV3",ARC3,IV3,SiO2);
				new G4LogicalBorderSurface("LDARC3",LD3,ARC3,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC3",fLMO1,ARC3,SiO);
				new G4LogicalBorderSurface("ARCLD3",ARC3,LD3,SiO);
				//new G4LogicalBorderSurface("ARCLMO3",ARC3,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC3",LD3,ARC3,SiO2);
			}
		} else if(face3 == 'e'){
			G4cout << "\nFace 3: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM3",fLMO1,EM3,Gs);
		}
	
		if(face4 == 'l'){ 
			G4cout << "\nFace 4: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC4",IV4,ARC4,SiO);
				new G4LogicalBorderSurface("ARCLD4",ARC4,LD4,SiO);
				//new G4LogicalBorderSurface("ARCIV4",ARC4,IV4,SiO2);
				new G4LogicalBorderSurface("LDARC4",LD4,ARC4,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC4",fLMO1,ARC4,SiO);
				new G4LogicalBorderSurface("ARCLD4",ARC4,LD4,SiO);
				//new G4LogicalBorderSurface("ARCLMO4",ARC4,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC4",LD4,ARC4,SiO2);
			}
		} else if(face4 == 'e'){
			G4cout << "\nFace 4: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM4",fLMO1,EM4,Gs);
		}
	
		if(face5 == 'l'){ 
			G4cout << "\nFace 5: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC5",IV5,ARC5,SiO);
				new G4LogicalBorderSurface("ARCLD5",ARC5,LD5,SiO);
				//new G4LogicalBorderSurface("ARCIV5",ARC5,IV5,SiO2);
				new G4LogicalBorderSurface("LDARC5",LD5,ARC5,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC5",fLMO1,ARC5,SiO);
				new G4LogicalBorderSurface("ARCLD5",ARC5,LD5,SiO);
				//new G4LogicalBorderSurface("ARCLMO5",ARC5,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC5",LD5,ARC5,SiO2);
			}
		} else if(face5 == 'e'){
			G4cout << "\nFace 5: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM5",fLMO1,EM5,Gs);
		}
	
		if(face6 == 'l'){ 
			G4cout << "\nFace 6: LD" << G4endl;
			if(LDgap > 0){
				//new G4LogicalBorderSurface("IVARC6",IV6,ARC6,SiO);
				new G4LogicalBorderSurface("ARCLD6",ARC6,LD6,SiO);
				//new G4LogicalBorderSurface("ARCIV6",ARC6,IV6,SiO2);
				new G4LogicalBorderSurface("LDARC6",LD6,ARC6,SiO2);
			} else {
				//new G4LogicalBorderSurface("LMOARC6",fLMO1,ARC6,SiO);
				new G4LogicalBorderSurface("ARCLD6",ARC6,LD6,SiO);
				//new G4LogicalBorderSurface("ARCLMO6",ARC6,fLMO1,SiO2);
				new G4LogicalBorderSurface("LDARC6",LD6,ARC6,SiO2);
			}
		} else if(face6 == 'e'){
			G4cout << "\nFace 6: EM" << G4endl;
		//	new G4LogicalBorderSurface("LMOEM6",fLMO1,EM6,Gs);
		}
	#endif
	}
#endif



	// *****************
	// MISC
	// *****************


#if ELD
	int dc = 10; // take into account largest gap between crys and ld/em

	G4double ELDPos[19]
		= {0, 0, LMOHalfSizes[2] + LDHalfSizes[2] + dc,
			0, 0, -LMOHalfSizes[2] - LDHalfSizes[2] - dc,
			0, -LMOHalfSizes[1] - LDHalfSizes[2] - dc, 0,
			0, LMOHalfSizes[1] + LDHalfSizes[2] + dc, 0,
			-LMOHalfSizes[0] - LDHalfSizes[2] - dc, 0, 0,
			LMOHalfSizes[0] + LDHalfSizes[2] + dc, 0, 0}; 

	// All LDHS' because rotated appropriately

	auto EncLD = new G4Box("EncLDs",ELDPos[2]*mm,ELDPos[2]*mm,0.5*mm);    

	G4LogicalVolume* EncLDs
		= new G4LogicalVolume(EncLD,LightDetectorMaterial,"EncLDs"); 


	//G4UserLimits* lim = new G4UserLimits();
	//lim->SetMaxAllowedStep(0.05*mm);
	//EncLDs->SetUserLimits(lim);


	G4VisAttributes* ELDVA = new G4VisAttributes();
	ELDVA->SetColor(0.0,0.5,0.5,0.5);
	ELDVA->SetForceSolid(true);
	EncLDs->SetVisAttributes(ELDVA);

	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);

	G4VPhysicalVolume* ELD1 = new G4PVPlacement(0,G4ThreeVector(ELDPos[0],ELDPos[1],ELDPos[2]),EncLDs,"ELD1",worldLV,false,80,fCheckOverlaps);	

	//G4VPhysicalVolume* ELD2 = new G4PVPlacement(0,G4ThreeVector(ELDPos[3],ELDPos[4],ELDPos[5]),EncLDs,"ELD2",worldLV,false,81,fCheckOverlaps);	

	G4VPhysicalVolume* ELD3 = new G4PVPlacement(t1,G4ThreeVector(ELDPos[6],ELDPos[7],ELDPos[8]),EncLDs,"ELD3",worldLV,false,82,fCheckOverlaps);	

	G4VPhysicalVolume* ELD4 = new G4PVPlacement(t1,G4ThreeVector(ELDPos[9],ELDPos[10],ELDPos[11]),EncLDs,"ELD4",worldLV,false,83,fCheckOverlaps);	

	G4VPhysicalVolume* ELD5 = new G4PVPlacement(p1,G4ThreeVector(ELDPos[12],ELDPos[13],ELDPos[14]),EncLDs,"ELD5",worldLV,false,84,fCheckOverlaps);	
	
	G4VPhysicalVolume* ELD6 = new G4PVPlacement(p1,G4ThreeVector(ELDPos[15],ELDPos[16],ELDPos[17]),EncLDs,"ELD6",worldLV,false,85,fCheckOverlaps);	



#endif

	//TODO: add back surfaces when particles are fixed






	// *****************
	// TEMP
	// *****************


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
	//TODO: I know I know cringe
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");
	char lmos = readExtFile("./options/lmos");
	char arc = readExtFile("./options/arcmat");

	G4VSensitiveDetector* fLDSD 
		= new LightDetectorSD("LightDetectorSD", 
				"AbsorberHitsCollection", fNofLayers);

	G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);

	G4VSensitiveDetector* fELDSD 
		= new ELightDetectorSD("ELightDetectorSD", 
				"HitsCollection", fNofLayers);
	
	G4SDManager::GetSDMpointer()->AddNewDetector(fELDSD);
	

	G4cout << "Crystal sensitive" << G4endl;
	SetSensitiveDetector("LMO_primary",fLDSD);

#if RegularGeometry
	if(face1 == 'l' || face2 == 'l' || face3 == 'l' || face4 == 'l' || face5 == 'l' || face6 == 'l'){	
	#if ARC
		if(arc == '2' || arc == '1' || arc == 'n'){
			SetSensitiveDetector("ARCLV",fLDSD); // LDs sensitive
			G4cout << "ARCs sensitive" << G4endl;
		}
	#endif
		SetSensitiveDetector("LDLV",fLDSD); // LDs sensitive
		G4cout << "LDs sensitive" << G4endl;
	}

	
	if(face1 == 'e' || face2 == 'e' || face3 == 'e' || face4 == 'e' || face5 == 'e' || face6 == 'e'){
	
		SetSensitiveDetector("EMLV",fLDSD);
		G4cout << "EMs sensitive" << G4endl;
	}


	// for troubleshooting missing abs peak
	//SetSensitiveDetector("World",fLDSD);
	
	#if ELD
	SetSensitiveDetector("EncLDs",fELDSD);
	#endif
#endif

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

