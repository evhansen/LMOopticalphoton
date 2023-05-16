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


/////////////////
// Defining some general materials (mostly from other systems)
/////////////////
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



/////////////////
// Defining the world
/////////////////
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
	/////////////////
	// First read options 
	/////////////////
	
	// Gap between the top of the frame and topmost LDs (default is 10mm)
	float tgap = stof(readExtFileLine("./options/TES_topgap"));
	
	// ARC toggle 
	char arc = readExtFile("./options/TES_arc");
	
	// Reflector toggle (underside of top slab)
	char refc = readExtFile("./options/TES_ref");
	
	// Cu can encapsulating system toggle
	char canc = readExtFile("./options/TES_can");
	
	// 70% Reflectivity toggle for Cu can
	char crefc = readExtFile("./options/TES_canref");



	/////////////////
	// Defining Materials for the physical volumes
	/////////////////

	G4Material* FrameMaterial;
	G4Material* LDMaterial;
	G4Material* ARCMaterial;
	G4Material* ReflectorMaterial;

	DefineTESMisc(&FrameMaterial, &LDMaterial, &ARCMaterial, &ReflectorMaterial);



	/////////////////
	// Defining lengths :)
	/////////////////

	G4double FRAMESIDE = (55.5/2)*mm;
	G4double SIDELONG = (135/2)*mm; // 127mm + 8mm top slab thick
	G4double TOPLONG = (102/2)*mm;
	G4double LDSIDE = (45/2)*mm;

	// Note: [0] is x-length, [1] is y-length, and [2] is z-length (thickness)
	G4double SideHalfSizes[3] = {SIDELONG + tgap, FRAMESIDE, 2*mm};
	G4double TopHalfSizes[3] = {TOPLONG, FRAMESIDE, 4*mm};
	G4double SlatHalfSizes[3] = {TOPLONG, FRAMESIDE, 3.5*mm};
	G4double LDHalfSizes[3] = {LDSIDE, LDSIDE, 0.25*mm}; //500 micrometres
	G4double ARCHalfSizes[3] = {LDSIDE, LDSIDE, (0.000069/2)*mm}; //69 nm
	G4double RefHalfSizes[3] = {TOPLONG, FRAMESIDE, 0.000005*mm}; //1 microm

	G4double CanThickHS = 1*mm;
	G4double CanSideBuffHS = 3*mm; //0.5*mm; // halfsize
       	G4double CanHeightBuffHS = 10*mm;
	G4double CanRadius = sqrt((TopHalfSizes[0]*TopHalfSizes[0]) + (TopHalfSizes[1]*TopHalfSizes[1])) + CanSideBuffHS;

	G4double CanSideHalfSizes[5] = {CanRadius, CanRadius + CanThickHS, SideHalfSizes[0]*mm + CanHeightBuffHS, 0.*rad, 2*pi*rad}; 
	G4double CanEndHalfSizes[5] = {0.*mm, CanSideHalfSizes[1], CanThickHS, CanSideHalfSizes[3], CanSideHalfSizes[4]};
	//inner radius, outer radius, height, start angle, spanning angle



	/////////////////
	// Defining solids
	/////////////////

	// Note:
	// G4Box : x, y, z
	// G4Tubs : inner radius, outer radius, height, 2xangles
	G4VSolid* SideS = new G4Box("SideS",
			SideHalfSizes[0],SideHalfSizes[1],SideHalfSizes[2]);
	G4VSolid* TopS = new G4Box("TopS",
			TopHalfSizes[0],TopHalfSizes[1],TopHalfSizes[2]);
	G4VSolid* LDS = new G4Box("LDS",
			LDHalfSizes[0],LDHalfSizes[1],LDHalfSizes[2]);

	G4VSolid* Slattemp = new G4Box("Slattemp",
			SlatHalfSizes[0],SlatHalfSizes[1],SlatHalfSizes[2]);
	G4VSolid* Slattemp2 = new G4SubtractionSolid("Slattemp2",
			Slattemp, LDS, 0, G4ThreeVector((2*mm)+LDHalfSizes[0],0,0));
	G4VSolid* SlatS = new G4SubtractionSolid("SlatS",
			Slattemp2, LDS, 0, G4ThreeVector(-(2*mm)-LDHalfSizes[0],0,0));

	G4VSolid* LDIVS = new G4Box("LDIVS",
			LDHalfSizes[0],LDHalfSizes[1],LDHalfSizes[2]);
	G4VSolid* ARCS = new G4Box("ARCS",
			ARCHalfSizes[0],ARCHalfSizes[1],ARCHalfSizes[2]);
	G4VSolid* ARCIVS = new G4Box("LDIVS",
			ARCHalfSizes[0],ARCHalfSizes[1],ARCHalfSizes[2]);
	
	G4VSolid* ReflectorS = new G4Box("ReflectorS",
			RefHalfSizes[0],RefHalfSizes[1],RefHalfSizes[2]);

	G4VSolid* CanSideS = new G4Tubs("CanSideS",CanSideHalfSizes[0],CanSideHalfSizes[1],
			CanSideHalfSizes[2],CanSideHalfSizes[3],CanSideHalfSizes[4]);
	G4VSolid* CanEndS = new G4Tubs("CanEndS",CanEndHalfSizes[0],CanEndHalfSizes[1],
			CanEndHalfSizes[2], CanEndHalfSizes[3], CanEndHalfSizes[4]);



	/////////////////
	// Defining logical volumes and colours
	/////////////////
	
	G4LogicalVolume* SideLV = new G4LogicalVolume(SideS,FrameMaterial,"SideLV");
	G4LogicalVolume* TopLV = new G4LogicalVolume(TopS,FrameMaterial,"TopLV");
	G4LogicalVolume* SlatLV = new G4LogicalVolume(SlatS,FrameMaterial,"SlatLV");
	G4LogicalVolume* LDLV = new G4LogicalVolume(LDS,LDMaterial,"LDLV");
	//G4LogicalVolume* LDIVLV = new G4LogicalVolume(LDIVS,WorldMaterial,"LDIVLV");
	G4LogicalVolume* ARCLV = new G4LogicalVolume(ARCS,ARCMaterial,"ARCLV");
	G4LogicalVolume* ReflectorLV = new G4LogicalVolume(ReflectorS,ReflectorMaterial,"ReflectorLV");
	G4LogicalVolume* CanSideLV = new G4LogicalVolume(CanSideS,FrameMaterial,"CanSideLV");
	G4LogicalVolume* CanEndLV = new G4LogicalVolume(CanEndS,FrameMaterial,"CanEndLV");


	G4VisAttributes* FrameVA = new G4VisAttributes();
	FrameVA->SetColor(1.0,0.0,1.0,0.10);
	FrameVA->SetForceSolid(true);

	G4VisAttributes* LDVA = new G4VisAttributes();
	LDVA->SetColor(0.0,1.0,1.0,0.15);
	LDVA->SetForceSolid(true);

	G4VisAttributes* ARCVA = new G4VisAttributes();
	ARCVA->SetColor(0.0,1.0,0.0,0.15);

	G4VisAttributes* ReflectorVA = new G4VisAttributes();
	ReflectorVA->SetColor(0.0,0.0,1.0,0.15);
	ReflectorVA->SetForceSolid(true);

	G4VisAttributes* CanVA = new G4VisAttributes();
	CanVA->SetColor(0.0,0.7,0.3,0.15);

	SideLV->SetVisAttributes(FrameVA);
	SlatLV->SetVisAttributes(FrameVA);
	TopLV->SetVisAttributes(FrameVA);
	LDLV->SetVisAttributes(LDVA);
	//LDIVLV->SetVisAttributes(IVVA);
	ARCLV->SetVisAttributes(ARCVA);
	ReflectorLV->SetVisAttributes(ReflectorVA);
	CanSideLV->SetVisAttributes(CanVA);
	CanEndLV->SetVisAttributes(CanVA);



	/////////////////
	// Defining physical volumes and related (like the rotation matrix :))
	/////////////////

	G4RotationMatrix* test = new G4RotationMatrix();
	test->rotateY(90.*deg);
	test->rotateX(180.*deg);



	// FRAME //

	//Commented out the pointers to the data since I use skins instead of borders :)
	//G4VPhysicalVolume* Top = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SideHalfSizes[0]*2 - TopHalfSizes[2]),TopLV,"TopPV",worldLV,false,0,fCheckOverlaps);

	//G4VPhysicalVolume* Side1 = 
		new G4PVPlacement(test,G4ThreeVector(TopHalfSizes[0], 0, SideHalfSizes[0]),SideLV,"Side1PV",worldLV,false,1,fCheckOverlaps);
	//G4VPhysicalVolume* Side2 = 
		new G4PVPlacement(test,G4ThreeVector(-TopHalfSizes[0], 0, SideHalfSizes[0]),SideLV,"Side2PV",worldLV,false,2,fCheckOverlaps);

	G4double SlatInit = 17*mm - SlatHalfSizes[2];
	G4double SlatHeightDiff = 25*mm;

	//G4VPhysicalVolume* Slat1 = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SlatInit),SlatLV,"Slat1PV",worldLV,false,3,fCheckOverlaps);
	//G4VPhysicalVolume* Slat2 = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SlatInit + SlatHeightDiff),SlatLV,"Slat2PV",worldLV,false,4,fCheckOverlaps);
	//G4VPhysicalVolume* Slat3 = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SlatInit + (SlatHeightDiff * 2)),SlatLV,"Slat3PV",worldLV,false,5,fCheckOverlaps);
	//G4VPhysicalVolume* Slat4 = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SlatInit + (SlatHeightDiff * 3)),SlatLV,"Slat4PV",worldLV,false,6,fCheckOverlaps);
	//G4VPhysicalVolume* Slat5 = 
		new G4PVPlacement(0,G4ThreeVector(0, 0, SlatInit + (SlatHeightDiff * 4)),SlatLV,"Slat5PV",worldLV,false,7,fCheckOverlaps);



	// LDs //
	
	G4double LDInit = SlatInit + SlatHalfSizes[2] + LDHalfSizes[2];
	G4double LDSide = 2*mm+LDHalfSizes[0];

	//G4VPhysicalVolume* LD1 = 
		new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit),LDLV,"LD1",worldLV,false,20,fCheckOverlaps);
	//G4VPhysicalVolume* LD2 = 
		new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit),LDLV,"LD2",worldLV,false,21,fCheckOverlaps);

	//G4VPhysicalVolume* LD3 = 
		new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + SlatHeightDiff),LDLV,"LD3",worldLV,false,22,fCheckOverlaps);
	//G4VPhysicalVolume* LD4 = 
		new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + SlatHeightDiff),LDLV,"LD4",worldLV,false,23,fCheckOverlaps);

	//G4VPhysicalVolume* LD5 = 
		new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 2)),LDLV,"LD5",worldLV,false,24,fCheckOverlaps);
	//G4VPhysicalVolume* LD6 = 
		new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 2)),LDLV,"LD6",worldLV,false,25,fCheckOverlaps);

	//G4VPhysicalVolume* LD7 = 
		new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 3)),LDLV,"LD7",worldLV,false,26,fCheckOverlaps);
	//G4VPhysicalVolume* LD8 = 
		new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 3)),LDLV,"LD8",worldLV,false,27,fCheckOverlaps);

	//G4VPhysicalVolume* LD9 = 
		new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 4)),LDLV,"LD9",worldLV,false,28,fCheckOverlaps);
	//G4VPhysicalVolume* LD10 = 
		new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 4)),LDLV,"LD10",worldLV,false,29,fCheckOverlaps);



	// ARCs //

	if(arc == '1') {
		G4double ARCbuf = ARCHalfSizes[2] + LDHalfSizes[2];

		//G4VPhysicalVolume* ARC1 = 
			new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit+ARCbuf),ARCLV,"ARC1",worldLV,false,30,fCheckOverlaps);	
		//G4VPhysicalVolume* ARC2 = 
			new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit+ARCbuf),ARCLV,"ARC2",worldLV,false,31,fCheckOverlaps);

		//G4VPhysicalVolume* ARC3 = 
			new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + SlatHeightDiff+ARCbuf),ARCLV,"ARC3",worldLV,false,32,fCheckOverlaps);
		//G4VPhysicalVolume* ARC4 = 
			new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + SlatHeightDiff+ARCbuf),ARCLV,"ARC4",worldLV,false,33,fCheckOverlaps);
	
		//G4VPhysicalVolume* ARC5 = 
			new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 2)+ARCbuf),ARCLV,"ARC5",worldLV,false,34,fCheckOverlaps);
		//G4VPhysicalVolume* ARC6 = 
			new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 2)+ARCbuf),ARCLV,"ARC6",worldLV,false,35,fCheckOverlaps);
	
		//G4VPhysicalVolume* ARC7 = 
			new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 3)+ARCbuf),ARCLV,"ARC7",worldLV,false,36,fCheckOverlaps);
		//G4VPhysicalVolume* ARC8 = 
			new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 3)+ARCbuf),ARCLV,"ARC8",worldLV,false,37,fCheckOverlaps);
	
		//G4VPhysicalVolume* ARC9 = 
			new G4PVPlacement(0,G4ThreeVector(LDSide,0,LDInit + (SlatHeightDiff * 4)+ARCbuf),ARCLV,"ARC9",worldLV,false,38,fCheckOverlaps);
		//G4VPhysicalVolume* ARC10 = 
			new G4PVPlacement(0,G4ThreeVector(-LDSide,0,LDInit + (SlatHeightDiff * 4)+ARCbuf),ARCLV,"ARC10",worldLV,false,39,fCheckOverlaps);

	} 

/* TODO: add back borders at some point if necessary!
	else {
		// ADD LDIV

		// LDIV 
		if(ivc == '1') {

		}
	}
*/



	// Reflectivities for the reflector and can //

	G4double energies[2] = {0.0001*eV, 99999999999*eV};
	G4double canref[2] = {0.9,0.9};
	G4double refref[2] = {0.7,0.7};



	// REFLECTOR //

	if(refc == '1') {
		G4MaterialPropertiesTable* RefMPT = new G4MaterialPropertiesTable();
		RefMPT->AddProperty("REFLECTIVITY",energies,refref,2);
		//RefMPT->AddProperty("TRANSMITTANCE",energies,reftrans,2);
	
		G4OpticalSurface* RefSurf = new G4OpticalSurface("ReflectorSurface");
		RefSurf->SetType(dielectric_metal);
		RefSurf->SetFinish(polished);
		//RefSurf->SetModel();
		RefSurf->SetMaterialPropertiesTable(RefMPT);

		//G4LogicalSkinSurface* RefSkin = 
			new G4LogicalSkinSurface("ReflectorSkin",ReflectorLV,RefSurf);

		//G4VPhysicalVolume* Reflector = 
			new G4PVPlacement(0,G4ThreeVector(0, 0, SideHalfSizes[0]*2 - 2*TopHalfSizes[2] - RefHalfSizes[2]),ReflectorLV,"Reflector",worldLV,false,10,fCheckOverlaps);

	}



	// Cu can //
	
	if(canc == '1') {

		// Cu can reflectivity //
		
		if(crefc == '1') {
			G4MaterialPropertiesTable* CanMPT = new G4MaterialPropertiesTable();
			CanMPT->AddProperty("REFLECTIVITY",energies,canref,2);
			//CanMPT->AddProperty("TRANSMITTANCE",energies,cantrans,2);
	
			G4OpticalSurface* CanSurf = new G4OpticalSurface("CanSurface");
			CanSurf->SetType(dielectric_metal);
			CanSurf->SetFinish(polished);
			//CanSurf->SetModel();
			CanSurf->SetMaterialPropertiesTable(CanMPT);
		
			//G4LogicalSkinSurface* CanSideSkin = 
				new G4LogicalSkinSurface("CanSideSkin",CanSideLV,CanSurf);
			//G4LogicalSkinSurface* CanEndSkin = 
				new G4LogicalSkinSurface("CanEndSkin",CanEndLV,CanSurf);
	
		}	

		//G4VPhysicalVolume* CanSide = 
			new G4PVPlacement(0,G4ThreeVector(0, 0, SideHalfSizes[0]),CanSideLV,"CanSide",worldLV,false,11,fCheckOverlaps);

		G4double CanEndbuf = 10*mm;

		//G4VPhysicalVolume* CanTop = 
			new G4PVPlacement(0,G4ThreeVector(0, 0, 2*SideHalfSizes[0] + CanEndHalfSizes[2] + CanEndbuf),CanEndLV,"CanTop",worldLV,false,12,fCheckOverlaps);
		//G4VPhysicalVolume* CanBottom = 
			new G4PVPlacement(0,G4ThreeVector(0, 0, -1*(CanEndHalfSizes[2] + CanEndbuf)),CanEndLV,"CanBottom",worldLV,false,13,fCheckOverlaps);

	}

}




void DetectorConstruction::ConstructSDandField()
{
	G4VSensitiveDetector* fLDSD 
		= new LightDetectorSD("LightDetectorSD", 
				"AbsorberHitsCollection", fNofLayers);

	G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);

	// No need for the sensitive detector that eats particles!
	//G4VSensitiveDetector* fELDSD 
	//	= new ELightDetectorSD("ELightDetectorSD", 
	//			"HitsCollection", fNofLayers);
	//G4SDManager::GetSDMpointer()->AddNewDetector(fELDSD);
	
	SetSensitiveDetector("LDLV",fLDSD);
}




/////////////////
// Read the first character from a file
/////////////////
char readExtFile(const char* filename){
	FILE *p = fopen(filename,"r");
	
	if(p == NULL){
		return '0';
	}
	
	char out = fgetc(p);

	fclose(p);

	return out;

}

/////////////////
// Read the first line from a file
/////////////////
std::string readExtFileLine(const char* filename){
	std::string out;
	std::ifstream p;

	p.open(filename);

	std::getline(p,out);

	return out;
}

