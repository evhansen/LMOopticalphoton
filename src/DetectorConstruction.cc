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
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

G4ThreadLocal

static char readExtFile(const char* filename);
std::string readExtFileLine(const char* filename);
float* readCSV(const char* filename, char delimiter,int*cols,int*rows);

	DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fCheckOverlaps(true),
fDetectorMessenger(nullptr),
fNofLayers(-1)
{
	// Set Options
	ldmat = readExtFile("./options/ldmat");
	emat = readExtFile("./options/emat");
	none = readExtFile("./options/none");
	sf = readExtFile("./options/surface");
	lmos = readExtFile("./options/lmos");
	oc = readExtFile("./options/oc");
	arc = readExtFile("./options/arc");
	pg = readExtFile("./options/pg");
	em = readExtFile("./options/em");
	dddm = readExtFile("./options/dddm");
	pLD5 = readExtFile("./options/pld5");


	ldg = stof(readExtFileLine("./options/ldg"));
	emg = ldg; //stof(readExtFileLine("./options/emg")); // TODO: make this different

	arct = stof(readExtFileLine("./options/arct"));
	
	fNofLayers = 1;
	
	EMLV = nullptr;

	DefineMaterials();
	
	fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{

	//delete fWorldMPT;
	delete worldLV;
	delete WorldMaterial;

	delete fLMOMPT;
	delete fLMO_LV_primary;
	delete LMOMaterial;

	delete fLDMPT;
	//delete LightDetectorLV;
	delete LightDetectorMaterial;

	//delete IVLV;
	
	delete EMLV;
	delete EMMaterial;

	//delete fSurfaceMPT;
	//delete fSurface;
	
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

	G4int numentries = 4;
	G4double energies[4] = {2.0*eV, 3.0*eV, 4.0*eV,500*keV};

	G4double a;  // mass of a mole;
	G4double z;  // z=mean number of protons;
	G4double density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;
	
	// Reading files
	int cols, rows, eit, nit, kit;
	float *vals;//, *e, *n, *k;


	// optical surface MPT (because of k -.-)
	LDk = new G4MaterialPropertiesTable();
	EMk = new G4MaterialPropertiesTable();
	
	// ********
	// World Material
	// ********
	
	new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);
	
	fWorldMPT   = new G4MaterialPropertiesTable();
	G4double rindex[4] = {1,1,1,1};
	fWorldMPT->AddProperty("RINDEX", energies, rindex, numentries);

	WorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Vacuum");
	WorldMaterial->SetMaterialPropertiesTable(fWorldMPT);

	if ( ! WorldMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nWORLDMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
	}

	// *******
	// LMO Materials
	// *******
	
	G4Element* elLi 
		= new G4Element(name="Lithium", symbol="Li", iz=3., a=6.94*g/mole);
	G4Element* elMo 
		= new G4Element(name="Molybdenum",symbol="Mo", iz=42, a=95.96*g/mole);
	G4Element* elO 
		= new G4Element(name="Oxygen", symbol="O", iz=8., a=16.00*g/mole);

	
	// 280 g, 45mm cube
	G4Material* Li2MoO4 =
		new G4Material (name="Li2MoO4", density=3.073*g/cm3, 3);
	
	Li2MoO4->AddElement(elLi, natoms=2);
	Li2MoO4->AddElement(elMo, natoms=1);
	Li2MoO4->AddElement(elO, natoms=4);

	G4double abslength[] = {100*mm, 100*mm, 100*mm, 100*mm};
	rindex[0] = 1.44;	rindex[1] = rindex[0];
	rindex[2] = rindex[0];	rindex[3] = rindex[0];

	fLMOMPT    = new G4MaterialPropertiesTable();	
	fLMOMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
	fLMOMPT->AddProperty("RINDEX", energies, rindex, numentries);

	LMOMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Li2MoO4");
	LMOMaterial->SetMaterialPropertiesTable(fLMOMPT);
	
	
	if (! LMOMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLMOMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
	}

	
	// ***********************
	// Energy-wavelength reference
	// ***********************
	
	// 2 eV ~ 0.6199 micrometres
	// 3 eV ~ 0.4133 micrometres
	// 4 eV ~ 0.31 micrometres
	// 500 keV ~ 2.4797*10^(-6) micrometers

	// **********
	// LD Materials
	// **********

	fLDMPT = new G4MaterialPropertiesTable();

	G4double trans[4] = {0.1,0.1,0.1,0.1};

	// apparently ... for dielectric-dielectric surfaces,
	// transmittance is for the proportion of g transmitted
	// through the material and the reflectivity is derived by
	// "R = 100% - T" (whatever that means). 

	G4double refl[4] = {0.9,0.9,0.9,0.9}; 

	// reflectivity specifies .. the absorption length since it ought be
	// R = 1 - a, where a is the fraction of g absorbed .. 
	
	// altogether, we have .1 transmitted, and .1 absorbed,
	// so that .8 undergoes refraction, reflection, TIR, etc.
	
	// for a dielectric_metal boundary,
	// reflectivity sets the probability for which g are reflected 
	// (exactly what we want so don't need to change)

	G4double Rrindex[4], Irindex[4];

	if(ldmat == '1'){
		new G4Material(name="Silicon", 
				z=14.0, a=28.0855*g/mole, density=2.33*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Silicon");

	
		vals = readCSV("./options/alrindex.csv",',',&cols,&rows);
		
		//e = (float*)malloc((rows+1)*sizeof(float*));
		//n = (float*)malloc((rows+1)*sizeof(float*));
		//k = (float*)malloc((rows+1)*sizeof(float*));
	
		eit = 0; G4double lden[rows] = {0};
		nit = 0; G4double ldrn[rows] = {0};
		kit = 0; G4double ldkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				//*(e+eit) = *(vals+i);
				//G4cout<<"Energy: "<<*(vals+i)<<G4endl;
				lden[eit] = 1239.8/((*(vals+i))*1000); 
				// en in micrometers
				eit++;
			}
	
			if((i+1) % 3 == 0){ // col 2
				//*(k+kit) = *(vals+i);
				//G4cout<<"Ext coeff k: "<<*(vals+i)<<G4endl;
				ldkn[kit] = *(vals+i);
				kit++;
			}
	
			if((i+2) % 3 == 0){ // col 1
				//*(n+nit) = *(vals+i);
				//G4cout<<"Rindx n: "<<*(vals+i)<<G4endl;
				ldrn[nit] = *(vals+i);
				nit++;
			}
		}
	


		//Last is actually at 0.2066 micrometers
		Rrindex[0] = 3.906;	Rrindex[1] = 5.222;
		Rrindex[2] = 5.01; 	Rrindex[3] = 1.01;	
		Irindex[0] = 0.022; 	Irindex[1] = 0.269;
		Irindex[2] = 3.586; 	Irindex[3] = 2.909;
		//rindex[0] = ;		rindex[1] = ;
		//rindex[2] = ;		rindex[3] = ;

		abslength[0] = pow(10,2)*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		/*	
		fLDMPT->AddProperty("RINDEX", 
				lden, ldrn, rows);
		fLDMPT->AddProperty("REALRINDEX", 
				lden, ldrn, rows);
		fLDMPT->AddProperty("IMAGINARYRINDEX", 
				lden, ldkn, rows);
		LDk->AddProperty("IMAGINARYRINDEX", 
				lden, ldkn, rows);
		*/

		fLDMPT->AddProperty("RINDEX", 
			energies, Rrindex, numentries);
		fLDMPT->AddProperty("REALRINDEX", 
			energies, Rrindex, numentries);
		fLDMPT->AddProperty("IMAGINARYRINDEX", 	
			energies, Irindex, numentries);
		LDk->AddProperty("IMAGINARYRINDEX", 
			energies,Irindex, numentries);
				


	} else {
	
		new G4Material(name="Germanium", 
				z=32.0, a=72.61*g/mole, density=5.323*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Germanium");


		vals = readCSV("./options/gerindex.csv",',',&cols,&rows);
	
		eit = 0; G4double lden[rows] = {0};
		nit = 0; G4double ldrn[rows] = {0};
		kit = 0; G4double ldkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				lden[eit] = 1239.8/((*(vals+i))*1000); 
				eit++;
			}
	
			if((i+1) % 3 == 0){ // col 2
				ldkn[kit] = *(vals+i);
				kit++;
			}
	
			if((i+2) % 3 == 0){ // col 1
				ldrn[nit] = *(vals+i);
				nit++;
			}
		}


		//Last is actually at 0.2066 micrometers
		//Rrindex[0] = 5.588; Rrindex[1] = 4.082;
		//Rrindex[2] = 3.905; Rrindex[3] = 1.023;
		//Irindex[0] = 0.933; Irindex[1] = 2.145;
		//Irindex[2] = 3.336; Irindex[3] = 2.774;
		//rindex[0] = ;		rindex[1] = ;
		//rindex[2] = ;		rindex[3] = ;
		
		abslength[0] = 10*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];
	
		fLDMPT->AddProperty("RINDEX", 
				lden, ldrn, rows);
		fLDMPT->AddProperty("REALRINDEX", 
				lden, ldrn, rows);
		fLDMPT->AddProperty("IMAGINARYRINDEX", 
				lden, ldkn, rows);
		LDk->AddProperty("IMAGINARYRINDEX", 
				lden, ldkn, rows);
	}

	fLDMPT->AddProperty("ABSLENGTH",
			energies, abslength, numentries);	
/*	fLDMPT->AddProperty("RINDEX", 
			energies, Rrindex, numentries);
	fLDMPT->AddProperty("REALRINDEX", 
			energies, Rrindex, numentries);
	fLDMPT->AddProperty("IMAGINARYRINDEX", 
			energies, Irindex, numentries);*/

	if(dddm == '1'){
		LDk->AddProperty("TRANSMITTANCE", 
			energies, trans, numentries);
		LDk->AddProperty("REFLECTIVITY", 
			energies, refl, numentries);
	} else if (dddm == '2') {
		LDk->AddProperty("REFLECTIVITY",
			energies,refl,numentries);
	}

	LightDetectorMaterial->SetMaterialPropertiesTable(fLDMPT);

	if (! LightDetectorMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLDMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
	}


	// ********************
	// Encapsulatory Material
	// ********************

	EMMPT = new G4MaterialPropertiesTable();
		
	if(emat == '1'){
	
		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Cu");

		vals = readCSV("./options/curindex.csv",',',&cols,&rows);
	
		eit = 0; G4double emen[rows] = {0};
		nit = 0; G4double emrn[rows] = {0};
		kit = 0; G4double emkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				emen[eit] = 1239.8/((*(vals+i))*1000); 
				eit++;
			}
	
			if((i+1) % 3 == 0){ // col 2
				emkn[kit] = *(vals+i);
				kit++;
			}
	
			if((i+2) % 3 == 0){ // col 1
				emrn[nit] = *(vals+i);
				nit++;
			}
		}


		// Last one is actually at 0.1879 micrometers		
		//Rrindex[0] = 0.29419; Rrindex[1] = 1.28;
		//Rrindex[2] =  1.3814; Rrindex[3] = 0.94;
		//Irindex[0] =  3.2443; Irindex[1] = 2.207;
		//Irindex[2] =  1.7254; Irindex[3] = 1.337;
		//rindex[0] = ;		rindex[1] = ;
		//rindex[2] = ;		rindex[3] = ;

		abslength[0] = 10*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];
	
		EMMPT->AddProperty("REALRINDEX",emen,emrn,rows);
		EMMPT->AddProperty("RINDEX",emen,emrn,rows);
		EMMPT->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
		EMk->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);

/*	}else if (emat == '2'){
	
		EMMaterial = G4NistManager::Instance()->FindOrBuildMaterial
			("G4_Pyrex_Glass",isotope);
	
	}else if (emat == '3'){
	
		EMMaterial = G4NistManager::Instance()->FindOrBuildMaterial
			("G4_GLASS_LEAD",isotope);
	
	} else if (emat == '4'){
	
		EMMaterial = G4NistManager::Instance()->FindOrBuildMaterial
			("G4_POLYTRIFLUOROCHLOROETHYLENE",isotope);
*/	
	} else {
	
		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Al");

		vals = readCSV("./options/alrindex.csv",',',&cols,&rows);
	
		eit = 0; G4double emen[rows] = {0};
		nit = 0; G4double emrn[rows] = {0};
		kit = 0; G4double emkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				emen[eit] = 1239.8/((*(vals+i))*1000); 
				eit++;
			}
	
			if((i+1) % 3 == 0){ // col 2
				emkn[kit] = *(vals+i);
				kit++;
			}
	
			if((i+2) % 3 == 0){ // col 1
				emrn[nit] = *(vals+i);
				nit++;
			}
		}



		// Last one is actually at 1.2399E-04 micrometers
		//Rrindex[0] = 1.3658; Rrindex[1] = 0.52141;
		//Rrindex[2] =0.28012; Rrindex[3] = 0.99999;
		//Irindex[0] = 7.4049; Irindex[1] =5.001;
		//Irindex[2] = 3.7086; Irindex[3] =8.241*pow(10,-8);

		//rindex[0] = ;		rindex[1] = ;
		//rindex[2] = ;		rindex[3] = ;
	
		abslength[0] = pow(10,3)*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];
	
		EMMPT->AddProperty("REALRINDEX",emen,emrn,rows);
		EMMPT->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
		EMMPT->AddProperty("RINDEX",emen,emrn,rows);
		EMk->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
	}

	EMMPT->AddProperty("ABSLENGTH",energies, abslength, numentries);
	/*EMMPT->AddProperty("REALRINDEX",energies,Rrindex,numentries);
	EMMPT->AddProperty("IMAGINARYRINDEX",energies,Irindex,numentries);
	EMMPT->AddProperty("RINDEX",energies,Rrindex,numentries);*/

	if(dddm == '1'){
		EMk->AddProperty("TRANSMITTANCE",energies,trans,numentries);
		EMk->AddProperty("REFLECTIVITY",energies,refl,numentries);
	} else if(dddm=='2') {
		EMk->AddProperty("REFLECTIVITY",energies,refl,numentries);
	}

	EMMaterial->SetMaterialPropertiesTable(EMMPT);

	if (!EMMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nEMMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
	}	

	
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;


	// *****************
	// Optical surface MPTs?
	// *****************
	
	// *****************
	// Antireflective Coating
	// *****************

	if(arc == '1'){
		
		G4Element* elSi = new G4Element
			(name="Silicon", symbol="Si", iz=14, a=28.0855*g/mole);

		G4Element* elO = new G4Element
			(name="Oxygen", symbol="O", iz=8., a=16.00*g/mole);

		G4Material* SiO2 =
			new G4Material (name="SiO2", density=2.196*g/cm3, 2);
	
		SiO2->AddElement(elSi, natoms=1);
		SiO2->AddElement(elO, natoms=2);

		abslength[0] = 1/(383.94)*cm;  abslength[1]= abslength[0];
		abslength[2] = abslength[0];	abslength[3] = abslength[0];

		rindex[0] = 1.44;	rindex[1] = rindex[0];
		rindex[2] = rindex[0];	rindex[3] = rindex[0];

		ARCMPT= new G4MaterialPropertiesTable();	
		ARCMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
		ARCMPT->AddProperty("RINDEX", energies, rindex, numentries);

		ARCMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("SiO2");
		ARCMaterial->SetMaterialPropertiesTable(fLMOMPT);

	}

}

G4VPhysicalVolume* DetectorConstruction::DefineWorld(){

	// ******
	// World
	// ******

	G4double worldSizeXY = 10. * m;

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

	
	// *******
	// Uhhhh something else ...
	// *******
	
	
	return worldPV;
}


void DetectorConstruction::DefineVolumes()
{	

	// ******
	// Define cube sizes,
	// {x, y, z}
	// ******
	
	G4double LMOHalfSizes[3]
		= {22.5*mm, 22.5*mm, 22.5*mm};

	fLMO_xy = LMOHalfSizes[0]; // what if [0] != [1]?

	G4double LDHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double EMHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], 0.25*mm};

	G4double ARCHalfSizes[3]
		= {LMOHalfSizes[0], LMOHalfSizes[1], arct*nm};

	G4double LDgap;

	if(oc == '1'){ // all in contact with crystal
	   	LDgap = 0.*mm;
	} else { //
		LDgap = ldg*mm;
	} 

	// TODO: make EM, LD gaps diff
	G4double EMgap = LDgap;


	// ********
	// Positions
	// ********

	G4double LDPos[19]
		= {0, 0, LMOHalfSizes[2] + LDgap + LDHalfSizes[2],
		0, 0, -LMOHalfSizes[2] - LDgap - LDHalfSizes[2],
		0, -LMOHalfSizes[1] - LDgap - LDHalfSizes[2], 0,
		0, LMOHalfSizes[1] + LDgap + LDHalfSizes[2], 0,
		-LMOHalfSizes[0] - LDgap - LDHalfSizes[2], 0, 0,
		LMOHalfSizes[0] + LDgap + LDHalfSizes[2], 0, 0};
	
	G4double EMPos[19]
		= {0, 0, LMOHalfSizes[2] + EMgap + EMHalfSizes[2],
		0, 0, -LMOHalfSizes[2] - EMgap - EMHalfSizes[2],
		0, -LMOHalfSizes[1] - EMgap - EMHalfSizes[2], 0,
		0, LMOHalfSizes[1] + EMgap + EMHalfSizes[2], 0,
		-LMOHalfSizes[0] - EMgap - EMHalfSizes[2], 0, 0,
		LMOHalfSizes[0] + EMgap + EMHalfSizes[2], 0, 0};

	if(arc == '1'){
		LDPos[2] += 2*ARCHalfSizes[2];
		LDPos[5] -= 2*ARCHalfSizes[2];
		LDPos[7] -= 2*ARCHalfSizes[2];
		LDPos[10] += 2*ARCHalfSizes[2];
		LDPos[12] -= 2*ARCHalfSizes[2];
		LDPos[15] += 2*ARCHalfSizes[2];
	}

	G4double ARCPos[19]
		= {0, 0, LMOHalfSizes[2] + LDgap + ARCHalfSizes[2],
		0, 0, -LMOHalfSizes[2] - LDgap - ARCHalfSizes[2],
		0, -LMOHalfSizes[1] - LDgap - ARCHalfSizes[2], 0,
		0, LMOHalfSizes[1] + LDgap + ARCHalfSizes[2], 0,
		-LMOHalfSizes[0] - LDgap - ARCHalfSizes[2], 0, 0,
		LMOHalfSizes[0] + LDgap + ARCHalfSizes[2], 0, 0};

	G4double IVHalfSizes[4]; // [2],[3] for LDgap and EMgap respectively
	
	// Want IV in contact with both LMO and other materials
	// so the halfsize should be the distance between 
	// IVpos and other material
		
	G4double IVPos[19] = {0};

	//IVPos = {0, 0, Pos[0] / 2,		What we want!
	// 0, 0, -Pos[0] / 2,		Need to figure out when
	// 0, -Pos[1] / 2, 0,		to use LD and EM pos.
	// 0, Pos[1] / 2, 0,		LMO is at 0,0,0
	// -Pos[2] / 2, 0, 0,
	// Pos[2] / 2, 0, 0};

	// note: if in optical contact, need to mark
	// gaps as nonzero in case we have 
	// no 1-4 LDs, EMs
	// and need the IV to def LMO border


	if (oc == '1') { // if in optical contact
		LDgap = 0.5*mm; EMgap = LDgap;
	}

	if((none == '1') || (em != '1')){ 
		// all LDs or (the last) two LDs only
	
		IVPos[2] = (LDgap / 2) + LMOHalfSizes[2]; 
		IVPos[5] = -(LDgap / 2) - LMOHalfSizes[2];
		IVPos[7] = -(LDgap / 2) - LMOHalfSizes[1]; 
		IVPos[10] = (LDgap / 2) + LMOHalfSizes[1];
		IVPos[12] = -(LDgap / 2) - LMOHalfSizes[0]; 
		IVPos[15] = (LDgap / 2) + LMOHalfSizes[0];
	
		IVHalfSizes[0] = LMOHalfSizes[0]; // x 
		IVHalfSizes[1] = LMOHalfSizes[1]; // y
	
		IVHalfSizes[2] = LDgap / 2; // z
		IVHalfSizes[3] = 0; // z

	} else { // two LDs and four EMs
		
		IVPos[2] = (EMgap / 2) + LMOHalfSizes[2]; 
		IVPos[5] = -(EMgap / 2) - LMOHalfSizes[2];
		IVPos[7] = -(EMgap / 2) - LMOHalfSizes[1]; 
		IVPos[10] = (EMgap / 2) + LMOHalfSizes[1];
		IVPos[12] = -(LDgap / 2) - LMOHalfSizes[0];
		IVPos[15] = (LDgap / 2) + LMOHalfSizes[0];
	
		IVHalfSizes[0] = LMOHalfSizes[0]; // x 
		IVHalfSizes[1] = LMOHalfSizes[1]; // y
		
		IVHalfSizes[2] = LDgap / 2; // z
		IVHalfSizes[3] = EMgap / 2; // z
	}

	G4double* pntOCP; // pointer to arrays


	// ********
	// Rotations
	// ********
	
	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);

		
	// *****************
	// Defining Volumes
	// *****************
	
	// ******
	// LMO
	// ******

	pntOCP = LMOHalfSizes; //Get LMO sizes

	auto LMO_box
	= new G4Box("LMO",          // its name
	*(pntOCP+0), *(pntOCP+1), *(pntOCP+2));    //its size
	
	fLMO_LV_primary
	= new G4LogicalVolume(
	LMO_box,         // its solid
	LMOMaterial,    // its material
	"LMO_primary");          // its name

	G4VisAttributes* LMOVA = new G4VisAttributes();
	LMOVA->SetColor(1.0,1.0,0.0,0.3);
	LMOVA->SetForceSolid(true);
	fLMO_LV_primary->SetVisAttributes(LMOVA);

	G4VPhysicalVolume* fLMO1;

	if(lmos == '1' || lmos == '2'){
		fLMO1 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(0,0,0),  // place off center eventually
		fLMO_LV_primary,          // its logical volume
		"LMO1",          // its name
		worldLV,                // its mother  volume
		false,            // no boolean operation
		1,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	// ******
	// Interim Vacuum
	// ******	
	
	G4VPhysicalVolume* IV1,*IV2,*IV3,*IV4,*IV5,*IV6;
	
	if(oc != '1') { // not in optical contact

		pntOCP = IVHalfSizes; //Get IV sizes

		// for IV between LMO and LD
		auto LDInterimVacuum    
			= new G4Box("LDInterimVacuum",     // its name
			*(pntOCP+0), *(pntOCP+1), *(pntOCP+2));

		G4LogicalVolume* LDIVLV
			= new G4LogicalVolume(
			LDInterimVacuum,     // its solid
			WorldMaterial,  // its material
		 	"LDInterimVacuumLV");   // its name

		G4VisAttributes* IVVA = new G4VisAttributes();
		IVVA->SetColor(1.0,0.0,1.0,0.15);
		IVVA->SetForceSolid(true);
		LDIVLV->SetVisAttributes(IVVA);
	
		
		if((none == '1') || (em != '1')){ // LDs, no EM

			pntOCP = IVPos; //Get IV position

			IV1 = new G4PVPlacement(
			 	0,                // no rotation
			 	G4ThreeVector(*(pntOCP+0),
					*(pntOCP+1),*(pntOCP+2)),
				LDIVLV,          // its logical volume
			 	"IV1T",    // its name
				worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	30,                // copy number
			 	fCheckOverlaps);  // checking overlaps

			IV2 = new G4PVPlacement(
			 	0,                // no rotation 
			 	G4ThreeVector(*(pntOCP+3),
					*(pntOCP+4),*(pntOCP+5)),
			 	LDIVLV,          // its logical volume
			 	"IV1B",    // its name
			 	worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	31,                // copy number
			 	fCheckOverlaps);  // checking overlaps

		      	IV3 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+6),
					*(pntOCP+7),*(pntOCP+8)),
			  	LDIVLV,          // its logical volume
			  	"IV2B",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	32,                // copy number
			  	fCheckOverlaps);  // checking overlaps
	
		      	IV4 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+9),
					*(pntOCP+10),*(pntOCP+11)),
			  	LDIVLV,          // its logical volume
			  	"IV3T",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	33,                // copy number
			  	fCheckOverlaps);  // checking overlaps

	
		 } else { // if there is EM

			pntOCP = IVHalfSizes; //Get IV sizes
			
			// for IV between LMO and EM
			auto EMInterimVacuum    
			= new G4Box("EMInterimVacuum",     // its name
			*(pntOCP+0), *(pntOCP+1), *(pntOCP+3));

			G4LogicalVolume* EMIVLV
			= new G4LogicalVolume(
			EMInterimVacuum,     // its solid
			WorldMaterial,  // its material
		 	"EMInterimVacuumLV");   // its name

			EMIVLV->SetVisAttributes(IVVA);

			pntOCP = IVPos; //Get IV position
			
			IV1 = new G4PVPlacement(
			 	0,                // no rotation
			 	G4ThreeVector(*(pntOCP+0),
					*(pntOCP+1),*(pntOCP+2)),
				EMIVLV,          // its logical volume
			 	"IV1T",    // its name
				worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	30,                // copy number
			 	fCheckOverlaps);  // checking overlaps

			IV2 = new G4PVPlacement(
			 	0,                // no rotation 
			 	G4ThreeVector(*(pntOCP+3),
					*(pntOCP+4),*(pntOCP+5)),
			 	EMIVLV,          // its logical volume
			 	"IV1B",    // its name
			 	worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	31,                // copy number
			 	fCheckOverlaps);  // checking overlaps

		      	IV3 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+6),
					*(pntOCP+7),*(pntOCP+8)),
			  	EMIVLV,          // its logical volume
			  	"IV2B",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	32,                // copy number
			  	fCheckOverlaps);  // checking overlaps
	
		      	IV4 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+9),
					*(pntOCP+10),*(pntOCP+11)),
			  	EMIVLV,          // its logical volume
			  	"IV3T",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	33,                // copy number
			  	fCheckOverlaps);  // checking overlaps
		} 

		//pntOCP = IVPos; //Get IV position
	    	
		// Will always have these two LDs	
		IV5 = new G4PVPlacement(
		  	p1,                // 
		 	G4ThreeVector(*(pntOCP+12),
				*(pntOCP+13),*(pntOCP+14)),
		  	LDIVLV,          // its logical volume
		  	"IV3B",    // its name
		  	worldLV,          // its mother  volume
		  	false,            // no boolean operation
		  	34,                // copy number
		  	fCheckOverlaps);  // checking overlaps

		IV6 = new G4PVPlacement(
		  	p1,                // 
		 	G4ThreeVector(*(pntOCP+15),
				*(pntOCP+16),*(pntOCP+17)),
		  	LDIVLV,          // its logical volume
		  	"IV3T",    // its name
		  	worldLV,          // its mother  volume
		  	false,            // no boolean operation
		  	35,                // copy number
		  	fCheckOverlaps);  // checking overlaps
	
	} else { // in optical contact
	
		if (none != '1' && em != '1') {
		// if not all LDs and no em, 
		// need IV to def LMO surface

			pntOCP = IVHalfSizes; //Get IV sizes

			// for IV def LMO surface
			auto LMOInterimVacuum    
				= new G4Box("LMOInterimVacuum",     // its name
				*(pntOCP+0), *(pntOCP+1), *(pntOCP+2));

			G4LogicalVolume* LMOIVLV
				= new G4LogicalVolume(
				LMOInterimVacuum,     // its solid
				WorldMaterial,  // its material
			 	"LMOInterimVacuumLV");   // its name

			G4VisAttributes* IVVA = new G4VisAttributes();
			IVVA->SetColor(1.0,0.0,1.0,0.15);
			IVVA->SetForceSolid(true);
			LMOIVLV->SetVisAttributes(IVVA);
				
			pntOCP = IVPos; //Get IV position

			IV1 = new G4PVPlacement(
			 	0,                // no rotation
			 	G4ThreeVector(*(pntOCP+0),
					*(pntOCP+1),*(pntOCP+2)),
				LMOIVLV,          // its logical volume
			 	"IV1T",    // its name
				worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	30,                // copy number
			 	fCheckOverlaps);  // checking overlaps

			IV2 = new G4PVPlacement(
			 	0,                // no rotation 
			 	G4ThreeVector(*(pntOCP+3),
					*(pntOCP+4),*(pntOCP+5)),
			 	LMOIVLV,          // its logical volume
			 	"IV1B",    // its name
			 	worldLV,          // its mother  volume
			 	false,            // no boolean operation
			 	31,                // copy number
			 	fCheckOverlaps);  // checking overlaps

			IV3 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+6),
					*(pntOCP+7),*(pntOCP+8)),
			  	LMOIVLV,          // its logical volume
			  	"IV2B",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	32,                // copy number
			  	fCheckOverlaps);  // checking overlaps

		      	IV4 = new G4PVPlacement(
			  	t1,                // 
			 	G4ThreeVector(*(pntOCP+9),
					*(pntOCP+10),*(pntOCP+11)),
			  	LMOIVLV,          // its logical volume
			  	"IV3T",    // its name
			  	worldLV,          // its mother  volume
			  	false,            // no boolean operation
			  	33,                // copy number
			  	fCheckOverlaps);  // checking overlaps
		}
	
	}


	// ******
	// LDs
	// ******

	pntOCP = LDHalfSizes; //Get LD sizes

	auto LightDetector
	= new G4Box("LightDetector",     // its name
	*(pntOCP+0), *(pntOCP+1), *(pntOCP+2)); // its size

	G4LogicalVolume* LightDetectorLV
	= new G4LogicalVolume(
	LightDetector,     // its solid
	LightDetectorMaterial,  // its material
	"LDLV");   // its name

	pntOCP = ARCHalfSizes;

	auto ARC = new G4Box("Antirefcoat",
		*(pntOCP+0),*(pntOCP+1),*(pntOCP+2));

	G4LogicalVolume* ARCLV
	= new G4LogicalVolume(
	ARC,     // its solid
	ARCMaterial,  // its material
	"ARCLV");   // its name

	G4VisAttributes* ARCVA = new G4VisAttributes();
	ARCVA->SetColor(1.0,1.0,1.0,0.5);
	ARCVA->SetForceSolid(true);
	ARCLV->SetVisAttributes(ARCVA);


	G4VPhysicalVolume* LD1,*LD2,*LD3,*LD4,*LD5,*LD6;
	G4VPhysicalVolume* ARC1, *ARC2, *ARC3, *ARC4, *ARC5, *ARC6;
	G4VPhysicalVolume* EM1,*EM2,*EM3,*EM4,*EM5;
	
	if(none == '1'){ 
	// no other encapsulating materials, all LDs

		pntOCP = LDPos; //Get LD positions

		G4VPhysicalVolume* LD1 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(*(pntOCP+0),
			*(pntOCP+1),*(pntOCP+2)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD1T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		10,                // copy number
		fCheckOverlaps);  // checking overlaps

		G4VPhysicalVolume* LD2 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(*(pntOCP+3),
			*(pntOCP+4),*(pntOCP+5)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD1B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		11,                // copy number
		fCheckOverlaps);  // checking overlaps

		G4VPhysicalVolume* LD3 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(*(pntOCP+6),
			*(pntOCP+7),*(pntOCP+8)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD2T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		12,                // copy number
		fCheckOverlaps);  // checking overlaps

		G4VPhysicalVolume* LD4 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(*(pntOCP+9),
			*(pntOCP+10),*(pntOCP+11)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD2B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		13,                // copy number
		fCheckOverlaps);  // checking overlaps
		
		if(arc == '1'){
			
			pntOCP = ARCPos;

			G4VPhysicalVolume* ARC1 = new G4PVPlacement(
				0,                // no rotation
				G4ThreeVector(*(pntOCP+0),
					*(pntOCP+1),*(pntOCP+2)),  // at (0,0,0)
				ARCLV,          // its logical volume
				"ARC1T",    // its name
				worldLV,          // its mother  volume
				false,            // no boolean operation
				40,                // copy number
				fCheckOverlaps);  // checking overlaps

			G4VPhysicalVolume* ARC2 = new G4PVPlacement(
				0,                // no rotation
				G4ThreeVector(*(pntOCP+3),
					*(pntOCP+4),*(pntOCP+5)),  // at (0,0,0)
				ARCLV,          // its logical volume
				"ARC1B",    // its name
				worldLV,          // its mother  volume
				false,            // no boolean operation
				41,                // copy number
				fCheckOverlaps);  // checking overlaps
		
			G4VPhysicalVolume* ARC3 = new G4PVPlacement(
				t1,                //
				G4ThreeVector(*(pntOCP+6),
					*(pntOCP+7),*(pntOCP+8)),  // at (0,0,0)
				ARCLV,          // its logical volume
				"ARC2T",    // its name
				worldLV,          // its mother  volume
				false,            // no boolean operation
				42,                // copy number
				fCheckOverlaps);  // checking overlaps
		
			G4VPhysicalVolume* ARC4 = new G4PVPlacement(
				t1,                //
				G4ThreeVector(*(pntOCP+9),
					*(pntOCP+10),*(pntOCP+11)),  // at (0,0,0)
				ARCLV,          // its logical volume
				"ARC2B",    // its name
				worldLV,          // its mother  volume
				false,            // no boolean operation
				43,                // copy number
				fCheckOverlaps);  // checking overlaps
		
		}


		
		if(pLD5 == '1'){
					
			pntOCP = EMHalfSizes; //Get EM sizes
			
			auto EM
			= new G4Box("EncapsulatingMaterial",
			*(pntOCP+0),*(pntOCP+1),*(pntOCP+2));
		
			EMLV = new G4LogicalVolume(
			EM,     // its solid
			EMMaterial,  // its material
			"EMLV");   // its name

			G4VisAttributes* EMVA = new G4VisAttributes();
			EMVA->SetColor(0.0,1.0,0.0,1.0);
			EMVA->SetForceLineSegmentsPerCircle(4);
			EMLV->SetVisAttributes(EMVA);
				
			pntOCP = EMPos; //Get EM positions

			G4VPhysicalVolume* EM5 = new G4PVPlacement(
			p1,                //
			G4ThreeVector(*(pntOCP+12),
				*(pntOCP+13),*(pntOCP+14)), 
			EMLV,          // its logical volume
			"EM3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			24,                // copy number
			fCheckOverlaps);  // checking overlaps

		} else {
		
			pntOCP = LDPos;
			
			G4VPhysicalVolume* LD5 = new G4PVPlacement(
			p1,                // no rotation
			G4ThreeVector(*(pntOCP+12),
				*(pntOCP+13),*(pntOCP+14)),  // at (0,0,0)
			LightDetectorLV,          // its logical volume
			"LD3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			14,                // copy number
			fCheckOverlaps);  // checking overlaps
	
			if(arc == '1'){

				pntOCP = ARCPos;
				
				G4VPhysicalVolume* ARC5 = new G4PVPlacement(
				p1,                // no rotation
				G4ThreeVector(*(pntOCP+12),
					*(pntOCP+13),*(pntOCP+14)),  // at (0,0,0)
				ARCLV,          // its logical volume
				"ARC3T",    // its name
				worldLV,          // its mother  volume
				false,            // no boolean operation
				44,                // copy number
				fCheckOverlaps);  // checking overlaps
			
			}

		}
		
		pntOCP = LDPos;

		G4VPhysicalVolume* LD6 = new G4PVPlacement(
		p1,                // no rotation
		G4ThreeVector(*(pntOCP+15),
			*(pntOCP+16),*(pntOCP+17)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD3B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		15,                // copy number
		fCheckOverlaps);  // checking overlaps

		if(arc=='1'){
			pntOCP = ARCPos;

			G4VPhysicalVolume* ARC6 = new G4PVPlacement(
			p1,                // no rotation
			G4ThreeVector(*(pntOCP+15),
				*(pntOCP+16),*(pntOCP+17)),  // at (0,0,0)
			ARCLV,          // its logical volume
			"ARC3B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			45,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
		
	} else {

		pntOCP = LDPos; //Get LD positions

		G4VPhysicalVolume* LD5 = new G4PVPlacement(
		p1,                // no rotation
		G4ThreeVector(*(pntOCP+12),
			*(pntOCP+13),*(pntOCP+14)),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD3T",    // its name
		worldLV,          // its mother  volume			
		false,            // no boolean operation
		14,                // copy number
		fCheckOverlaps);  // checking overlaps
		
		if(arc == '1'){
			pntOCP = ARCPos; //Get LD positions

			G4VPhysicalVolume* ARC5 = new G4PVPlacement(
			p1,                // no rotation
			G4ThreeVector(*(pntOCP+12),
				*(pntOCP+13),*(pntOCP+14)),  // at (0,0,0)
			ARCLV,          // its logical volume
			"ARC3T",    // its name
			worldLV,          // its mother  volume			
			false,            // no boolean operation
			44,                // copy number
			fCheckOverlaps);  // checking overlaps
		}

	
		if(pg != '1'){	

			pntOCP = LDPos; //Get LD positions
			
			G4VPhysicalVolume* LD6 = new G4PVPlacement(
			p1,                // no rotation
			G4ThreeVector(*(pntOCP+15),
				*(pntOCP+16),*(pntOCP+17)), // at (0,0,0)
			LightDetectorLV,          // its logical volume
			"LD3B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			15,                // copy number
			fCheckOverlaps);  // checking overlaps

			if(arc == '1'){

				pntOCP = ARCPos; //Get LD positions
			
				G4VPhysicalVolume* ARC6 = new G4PVPlacement(
					p1,                // no rotation
					G4ThreeVector(*(pntOCP+15),
						*(pntOCP+16),*(pntOCP+17)),
					ARCLV,     
					"ARC3B",    // its name
					worldLV,          // its mother  volume
					false,            // no boolean operation
					45,                // copy number
					fCheckOverlaps);  // checking overlaps
			}


		}

		// ******
		// Photon gun hole
		// ******

		//if(false){
		if(pg == '1') { // using photon gun

			
			pntOCP = LDHalfSizes;

			auto LDHoles
			= new G4Tubs("LDHoles",
				0.0,0.25*mm,*(pntOCP+2),
				0.0,CLHEP::twopi); 

			G4SubtractionSolid* LDCH
			= new G4SubtractionSolid("LDCHoles",
				LightDetector,LDHoles);

			G4LogicalVolume* LDHLV
			= new G4LogicalVolume(
			LDCH,     // its solid
			WorldMaterial,  // its material
			"LDHolestLV");   // its name
			
			pntOCP = LDPos;

			G4VPhysicalVolume * LDHPV
			= new G4PVPlacement(
			p1,                // no rotation
			G4ThreeVector(*(pntOCP+15),
				*(pntOCP+16),*(pntOCP+17)),
			//G4ThreeVector(0,0,0),
			LDHLV,          // its logical volume
			"LDH1M",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			50,                // copy number
			fCheckOverlaps);  // checking overlaps
		
			G4Tubs* OFCyl = new G4Tubs("OF", 0.26*mm, 0.5*mm,
					20*mm, 0.*deg, 360.*deg);

			G4LogicalVolume* OFLV
			= new G4LogicalVolume(
			OFCyl,     // its solid
			WorldMaterial,  // its material
			"OFLV");   // its name

			G4VisAttributes* OFVA = new G4VisAttributes();
			OFVA->SetColor(0.2,0.5,1.0,0.5);
			OFVA->SetForceSolid(true);
			OFLV->SetVisAttributes(OFVA);


			G4VPhysicalVolume* OF 
			= new G4PVPlacement(p1,
			G4ThreeVector(22.5*mm+20*mm+2.5*mm,0,0),
			OFLV, "OFH1M", worldLV, false, 60, fCheckOverlaps);

			if(arc == '1'){

				pntOCP = ARCPos; //Get LD positions
			
				G4VPhysicalVolume* ARC6 = new G4PVPlacement(
					p1,                // no rotation
					G4ThreeVector(*(pntOCP+15),
						*(pntOCP+16),*(pntOCP+17)),
					ARCLV,     
					"ARC3B",    // its name
					worldLV,          // its mother  volume
					false,            // no boolean operation
					45,                // copy number
					fCheckOverlaps);  // checking overlaps
			}
		
		}


		// ******
		// Encapsulatory Material
		// ******
		
		if(em == '1'){
			
			pntOCP = EMHalfSizes; //Get EM sizes
			
			auto EM
			= new G4Box("EncapsulatingMaterial",
			*(pntOCP+0),*(pntOCP+1),*(pntOCP+2));
		
			EMLV = new G4LogicalVolume(
			EM,     // its solid
			EMMaterial,  // its material
			"EMLV");   // its name

			G4VisAttributes* EMVA = new G4VisAttributes();
			EMVA->SetColor(0.0,1.0,0.0,1.0);
			EMVA->SetForceLineSegmentsPerCircle(4);
			EMLV->SetVisAttributes(EMVA);
	
			pntOCP = EMPos; //Get EM positions

			
			G4VPhysicalVolume* EM1 = new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(*(pntOCP+0),
				*(pntOCP+1),*(pntOCP+2)),
			EMLV,          
			"EM1T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			20,                // copy number
			fCheckOverlaps);  // checking overlaps

			G4VPhysicalVolume* EM2 = new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(*(pntOCP+3),
				*(pntOCP+4),*(pntOCP+5)),  
			EMLV,          // its logical volume
			"EM1B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			21,                // copy number
			fCheckOverlaps);  // checking overlaps

			G4VPhysicalVolume* EM3 = new G4PVPlacement(
			t1,                //
			G4ThreeVector(*(pntOCP+6),
				*(pntOCP+7),*(pntOCP+8)),
			EMLV,          // its logical volume
			"EM2T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			22,                // copy number
			fCheckOverlaps);  // checking overlaps

			G4VPhysicalVolume* EM4 = new G4PVPlacement(
			t1,                //
			G4ThreeVector(*(pntOCP+9),
				*(pntOCP+10),*(pntOCP+11)), 
			EMLV,          // its logical volume
			"EM2B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			23,                // copy number
			fCheckOverlaps);  // checking overlaps

		}
	}


	// *****************
	// Defining Surfaces
	// *****************

	
	// have LDk and EMk
	G4OpticalSurface* ddg = new G4OpticalSurface("ddg");
	ddg->SetModel(unified);
	ddg->SetType(dielectric_dielectric);
	ddg->SetFinish(ground);
	ddg->SetSigmaAlpha(0.1);

	G4OpticalSurface* ddgEM = new G4OpticalSurface("ddgEM");
	ddgEM->SetModel(unified);
	ddgEM->SetType(dielectric_dielectric);
	ddgEM->SetFinish(ground);
	ddgEM->SetSigmaAlpha(0.1);
	ddgEM->SetMaterialPropertiesTable(EMk);

	G4OpticalSurface* ddgLD = new G4OpticalSurface("ddgLD");
	ddgLD->SetModel(unified);
	ddgLD->SetType(dielectric_dielectric);
	ddgLD->SetFinish(ground);
	ddgLD->SetSigmaAlpha(0.1);
	ddgLD->SetMaterialPropertiesTable(LDk);

	G4OpticalSurface* ddp = new G4OpticalSurface("ddp");
	ddp->SetModel(unified);
	ddp->SetType(dielectric_dielectric);
	ddp->SetFinish(polished);

	G4OpticalSurface* ddpEM = new G4OpticalSurface("ddpEM");
	ddpEM->SetModel(unified);
	ddpEM->SetType(dielectric_dielectric);
	ddpEM->SetFinish(polished);
	ddpEM->SetMaterialPropertiesTable(EMk);
	
	G4OpticalSurface* ddpLD = new G4OpticalSurface("ddpLD");
	ddpLD->SetModel(unified);
	ddpLD->SetType(dielectric_dielectric);
	ddpLD->SetFinish(polished);
	ddpLD->SetMaterialPropertiesTable(LDk);
	
	G4OpticalSurface* dmgEM = new G4OpticalSurface("dmgEM");
	dmgEM->SetModel(glisur);
	dmgEM->SetType(dielectric_metal);
	dmgEM->SetFinish(polished);
	dmgEM->SetMaterialPropertiesTable(EMk);

	G4OpticalSurface* dmgLD = new G4OpticalSurface("dmgLD");
	dmgLD->SetModel(glisur);
	dmgLD->SetType(dielectric_metal);
	dmgLD->SetFinish(polished);
	dmgLD->SetMaterialPropertiesTable(LDk);

	G4OpticalSurface* dmpEM = new G4OpticalSurface("dmpEM");
	dmpEM->SetModel(glisur);
	dmpEM->SetType(dielectric_metal);
	dmpEM->SetFinish(polished);
	dmpEM->SetMaterialPropertiesTable(EMk);
	
	G4OpticalSurface* dmpLD = new G4OpticalSurface("dmpLD");
	dmpLD->SetModel(glisur);
	dmpLD->SetType(dielectric_metal);
	dmpLD->SetFinish(polished);
	dmpLD->SetMaterialPropertiesTable(LDk);

	//dynamic_cast<G4OpticalSurface*>(surface1->GetSurface(fLMO1, worldPV)->GetSurfaceProperty());


	G4OpticalSurface *LMOIV, *IVLMO, *EMIV, *IVEM, *LDIV, *IVLD; 
	G4OpticalSurface *LMOEM, *EMLMO, *LMOLD, *LDLMO;

	/*if(sf == '3'){ // POLISHED
		sides = t3; tb = t3;
		G4cout << "POLISHED" << G4endl;

	} else if (sf == '2'){ // GROUND
		sides = t2; tb = t2;
		G4cout << "GROUND" << G4endl;

	} else { // BOTH */

		if(dddm=='1'){
			IVLD = ddpLD;	IVEM = ddpEM;
			LMOLD = ddgLD;	LMOEM = ddgEM;
		} else if (dddm=='2') {
			IVLD = dmpLD;	IVEM = dmpEM;
			LMOLD = dmgLD;	LMOEM = dmgEM;
		} else {
			IVLD = ddpLD;	IVEM = ddpEM;
			LMOLD = ddgLD;	LMOEM = ddgEM;
		}

		// polished
		LDIV = ddp; EMIV = ddp;
		EMLMO = ddp; LDLMO = ddp;

		// ground
		LMOIV = ddg; IVLMO = ddg; 

		G4cout << "BOTH" << G4endl;
	//}


	if(oc != '1'){ 
	// If there's an interim vacuum (not in optical contact)

	// TODO: check if LDgap and EMgap differ, i.e. optical contact diff

		if(lmos == '1' || lmos == '2'){
			new G4LogicalBorderSurface("LMOIV1",
				fLMO1,IV1,LMOIV);	
			new G4LogicalBorderSurface("LMOIV2",
				fLMO1,IV2,LMOIV);
			new G4LogicalBorderSurface("LMOIV3",
				fLMO1,IV3,LMOIV);
			new G4LogicalBorderSurface("LMOIV4",
				fLMO1,IV4,LMOIV);
			new G4LogicalBorderSurface("LMOIV5",
				fLMO1,IV5,LMOIV);
			new G4LogicalBorderSurface("LMOIV6",
				fLMO1,IV6,LMOIV);

			new G4LogicalBorderSurface("IVLMO1",
				IV1,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO2",
				IV2,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO3",
				IV3,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO4",
				IV4,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO5",
				IV5,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO6",
				IV6,fLMO1,IVLMO);
		}

		// always have two LDs
		new G4LogicalBorderSurface("IVLD6",
			IV6,LD6,IVLD);

		new G4LogicalBorderSurface("LDIV6",
			LD6,IV6,LDIV);
		
		if(pLD5 == '1'){
			new G4LogicalBorderSurface("IVEM5",
				IV5,EM5,IVEM);
	
			new G4LogicalBorderSurface("EMIVM5",
				EM5,IV5,EMIV);

		} else {
			new G4LogicalBorderSurface("IVLD5",
				IV5,LD5,IVLD);
	
			new G4LogicalBorderSurface("LDIV5",
				LD5,IV5,LDIV);
		}


		if(none == '1'){ // only LDs
			new G4LogicalBorderSurface("IVLD1",
				IV1,LD1,IVLD);
			new G4LogicalBorderSurface("IVLD2",
				IV2,LD2,IVLD);
			new G4LogicalBorderSurface("IVLD3",
				IV3,LD3,IVLD);
			new G4LogicalBorderSurface("IVLD4",
				IV4,LD4,IVLD);

			new G4LogicalBorderSurface("LDIV1",
				LD1,IV1,LDIV);
			new G4LogicalBorderSurface("LDIV2",
				LD2,IV2,LDIV);
			new G4LogicalBorderSurface("LDIV3",
				LD3,IV3,LDIV);
			new G4LogicalBorderSurface("LDIV4",
				LD4,IV4,LDIV);

		}else if(em == '1'){ // if there're EM
			new G4LogicalBorderSurface("EMIV1",
				EM1,IV1,IVLD);
			new G4LogicalBorderSurface("EMIV2",
				EM2,IV2,IVLD);
			new G4LogicalBorderSurface("EMIV3",
				EM3,IV3,IVLD);
			new G4LogicalBorderSurface("EMIV4",
				EM4,IV4,IVLD);

			new G4LogicalBorderSurface("IVEM1",
				IV1,EM1,LDIV);
			new G4LogicalBorderSurface("IVEM2",
				IV2,EM2,LDIV);
			new G4LogicalBorderSurface("IVEM3",
				IV3,EM3,LDIV);
			new G4LogicalBorderSurface("IVEM4",
				IV4,EM4,LDIV);
		} // else, neither! IV already def'd on LMO so good

	}else if(lmos == 1 || lmos == '2'){ 
	// If there's optical contact

		// always have two LDs
		new G4LogicalBorderSurface("LMOLD6",
			fLMO1,LD6,LMOLD);

		new G4LogicalBorderSurface("LDLMO6",
			LD6,fLMO1,LDLMO);

		
		if(pLD5 == '1'){
			new G4LogicalBorderSurface("LMOEM5",
				fLMO1,EM5,IVEM);
	
			new G4LogicalBorderSurface("EMLMO5",
				EM5,fLMO1,EMIV);

		} else {
		
			new G4LogicalBorderSurface("LMOLD5",
				fLMO1,LD5,LMOLD);
		
			new G4LogicalBorderSurface("LDLMO5",
				LD5,fLMO1,LDLMO);
		}
		
		
		
		if(none == '1'){ // all LDs
		
			new G4LogicalBorderSurface("LMOLD1",
				fLMO1,LD1,LMOLD);
			new G4LogicalBorderSurface("LMOLD2",
				fLMO1,LD2,LMOLD);
			new G4LogicalBorderSurface("LMOLD3",
				fLMO1,LD3,LMOLD);
			new G4LogicalBorderSurface("LMOLD4",
				fLMO1,LD4,LMOLD);

			new G4LogicalBorderSurface("LDLMO1",
				LD1,fLMO1,LDLMO);
			new G4LogicalBorderSurface("LDLMO2",
				LD2,fLMO1,LDLMO);
			new G4LogicalBorderSurface("LDLMO3",
				LD3,fLMO1,LDLMO);
			new G4LogicalBorderSurface("LDLMO4",
				LD4,fLMO1,LDLMO);
			
		} else if (em == '1') { // with four EMs

			new G4LogicalBorderSurface("LMOEM1",
				fLMO1,EM1,LMOEM);
			new G4LogicalBorderSurface("LMOEM2",
				fLMO1,EM2,LMOEM);
			new G4LogicalBorderSurface("LMOEM3",
				fLMO1,EM3,LMOEM);
			new G4LogicalBorderSurface("LMOEM4",
				fLMO1,EM4,LMOEM);

			new G4LogicalBorderSurface("EMLMO1",
				EM1,fLMO1,EMLMO);
			new G4LogicalBorderSurface("EMLMO2",
				EM2,fLMO1,EMLMO);
			new G4LogicalBorderSurface("EMLMO3",
				EM3,fLMO1,EMLMO);
			new G4LogicalBorderSurface("EMLMO4",
				EM4,fLMO1,EMLMO);
		
		}else { 
		// else, neither LDs or EM but still need IV to def LMO
			new G4LogicalBorderSurface("LMOIV1",
				fLMO1,IV1,LMOIV);	
			new G4LogicalBorderSurface("LMOIV2",
				fLMO1,IV2,LMOIV);
			new G4LogicalBorderSurface("LMOIV3",
				fLMO1,IV3,LMOIV);
			new G4LogicalBorderSurface("LMOIV4",
				fLMO1,IV4,LMOIV);

			new G4LogicalBorderSurface("IVLMO1",
				IV1,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO2",
				IV2,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO3",
				IV3,fLMO1,IVLMO);
			new G4LogicalBorderSurface("IVLMO4",
				IV4,fLMO1,IVLMO);
		}
	}



	// ********
	// Not sure
	// ********

	G4cout << G4endl << "Setting Optical Surface at Junctions." << G4endl;

	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
	auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	simpleBoxVisAtt->SetVisibility(true);
	LightDetectorLV->SetVisAttributes(simpleBoxVisAtt);

}


void DetectorConstruction::ConstructSDandField()
{
	//auto fLDSD
	G4VSensitiveDetector* fLDSD 
		= new LightDetectorSD("LightDetectorSD", 
				"AbsorberHitsCollection", fNofLayers);

	G4SDManager::GetSDMpointer()->AddNewDetector(fLDSD);


	SetSensitiveDetector("LDLV",fLDSD); // LDs sensitive
	
	if(false){ // exceptions when LV not init'd >:(
	//if(em == '1'){
		SetSensitiveDetector("EMLV",fLDSD);
		G4cout << "Encapsulatory Material sensitive" << G4endl;
	} else{
		G4cout << "Encapsulatory Material not sensitive" << G4endl;
	}


	if(true){
	//if(lmos == '2'){
		SetSensitiveDetector("LMO_primary",fLDSD);
		G4cout << "LMO sensitive" << G4endl;
	} else {
	
		G4cout << "LMO not sensitive" << G4endl;
	}
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


//int* readCSV(const char* filename, const char* delimiter){
float* readCSV(const char* filename, char delimiter,int*cols,int*rows){
	FILE *p = fopen(filename,"r");
	char *split = NULL;

	if(p == NULL){
		//return "0";
		G4cout << "File unreachable" << G4endl;
		return NULL;
	}

	char te;
	(*rows) = 0;
	(*cols) = 0;

	while((te = fgetc(p)) != EOF){
		if((*rows) == 0){
			if(te == delimiter){
				(*cols)++; 
				// x, y, z for [0], [1], [2]
				// so perfect correspondence
			}
		}

		if(te == '\n'){
			(*rows)++;
		}
	}

	//G4cout << "Rows: " << (*rows) << "\nCols: " << (*cols) << G4endl;

	rewind(p);

	float *arrpnts = (float*)malloc((*rows)*(*cols+1)*sizeof(float*));

	char out[200]; 
	int crows = 0;
	int ccols = 0;

	fgets(out,200,p);
	char * out_s = out + 3; // strip BOM

	do {
		out_s[strlen(out_s)-1] = '\0'; // really can't have line > 200
		split = strtok(out_s,&delimiter);

		while(split){

			if(ccols > (*cols)){
				break;	
			}

			//G4cout << split << " : " << ((((*cols)+1)*crows)+ccols) << G4endl;
			*(arrpnts + ((((*cols)+1)*crows)+ccols)) = atof(split);
			ccols++;
		
			split = strtok(NULL,&delimiter);
		}
	
		ccols = 0;
		crows++;
	
	} while(fgets(out_s,200,p));

	//for(int i = 0; i < ((*cols)+1)*(*rows); i++){
	//	G4cout << "Vals " << i << ": " << *(arrpnts+i) << G4endl;
	//}
	
	fclose(p);
	//free(arrpnts);

	return arrpnts;
}
