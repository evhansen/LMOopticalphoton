

/////////////////////////////////////////
// This is for defining materials 
// as one uses in DetectorConstruction
/////////////////////////////////////////


// File BOM issues
static int BOM = 0; // 0 if BOM needs to be removed, else do not remove.

// Read first char in file //
static char readExtFile(const char* filename);

// Read CSV //
float* readCSV(const char* filename, char delimiter,int*cols,int*rows);

// Set varying energy property using CSV values from readCSV //
void SetEnergyVaryingProperty(const char* filename, const char* property, const char delimiter, double energyunit, double propertyunit, G4MaterialPropertiesTable* MPT);

// Set varying energy properties using CSV values from readCSV //
void SetEnergyVaryingProperties(const char* filename, const char* property1, const char* property2, const char delimiter, double energyunit, double propertyunit1, double propertyunit2, G4MaterialPropertiesTable* MPT);



/////////////////
// World material and geometry
/////////////////

void DetectorConstruction::DefineWorldM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {1.0*eV, 3.0*eV, 4.0*eV,100*GeV};
	G4double a,z,density;
	G4String name= "";
	G4String symbol = "";


	new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);

	G4MaterialPropertiesTable* fWorldMPT = new G4MaterialPropertiesTable();

	G4double rindex[4] = {1,1,1,1}; 
	G4double abslength[4]={(200)*mm,(200)*mm,(200)*mm,(200)*mm};

	fWorldMPT->AddProperty("RINDEX", energies, rindex, numentries);
	fWorldMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
	
	WorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Vacuum");
	WorldMaterial->SetMaterialPropertiesTable(fWorldMPT);


	
	if ( ! WorldMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nWORLDMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}

	G4cout << "PRINTING WORLD TABLE\n" << G4endl;
	fWorldMPT->DumpTable();
}



/////////////////
// Crystal material (LMO and CoCalib runs)
/////////////////

void DetectorConstruction::DefineCrysM(char mat){
	G4double a,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	G4MaterialPropertiesTable* fLMOMPT = new G4MaterialPropertiesTable();	

	G4cout << "\nCrystal material: " << mat << "\n" << G4endl;

	if(mat == 't'){

		G4Element* elTe = G4NistManager::Instance()
			->FindOrBuildElement(52, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		// tetragonal	
		G4Material* TeO2 = new G4Material 
			(name="TeO2", density=6.04*g/cm3, 2);

		TeO2->AddElement(elTe, natoms=1);
		TeO2->AddElement(elO, natoms=2);

		LMOMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("TeO2");

		BOM = 2;

		SetEnergyVaryingProperty("./options/mpt/teo2-abs.csv", "ABSLENGTH", ',', 
				eV, cm, fLMOMPT);
	
		if(readExtFile("./options/teo2rndx") == 'o'){
			SetEnergyVaryingProperty("./options/mpt/teo2-norindex.csv", "RINDEX", ',', 
				eV, 1, fLMOMPT);
		} else {
			SetEnergyVaryingProperty("./options/mpt/teo2-nerindex.csv", "RINDEX", ',', 
				eV, 1, fLMOMPT);
		}	
	
		BOM = 0;

	} else if (mat == 'c') {

		/*
		G4Element* elCd = G4NistManager::Instance()
			->FindOrBuildElement(48, false);
		G4Element* elW = G4NistManager::Instance()
			->FindOrBuildElement(74, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		// tetragonal	
		G4Material* CdWO4 = new G4Material 
			(name="CdWO4", density=7.9*g/cm3, 6);

		CdWO4->AddElement(elCd, natoms=1);
		CdWO4->AddElement(elW, natoms=1);
		CdWO4->AddElement(elO, natoms=4);

		LMOMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("CdWO4");
		*/
	

		LMOMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_CdWO4");

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/cdwo4-rindex.csv", "RINDEX", ',', 
				eV, 1, fLMOMPT);
		SetEnergyVaryingProperty("./options/mpt/cdwo4-abs.csv", "ABSLENGTH", ',', 
				eV, cm, fLMOMPT);
		BOM = 0;

		std::vector<G4double> ener = {2.6*eV};
		std::vector<G4double> scin = {1200.};

		fLMOMPT->AddProperty("SCINTILLATIONCOMPONENT1",ener,scin);
		fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.005*ns);

		fLMOMPT->AddConstProperty("SCINTILLATIONYIELD",1200./keV);
		
		fLMOMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
		//fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.*ns);

	} else if (mat == 'p') {

		/*
		G4Element* elPb = G4NistManager::Instance()
			->FindOrBuildElement(82, false);
		G4Element* elW = G4NistManager::Instance()
			->FindOrBuildElement(74, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		// tetragonal	
		G4Material* PbWO4 = new G4Material 
			(name="PbWO4", density=8.28*g/cm3, 6);

		PbWO4->AddElement(elPb, natoms=1);
		PbWO4->AddElement(elW, natoms=1);
		PbWO4->AddElement(elO, natoms=4);

		LMOMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("PbWO4");
		*/

		LMOMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_PbWO4");

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/pbwo4-rindex.csv", "RINDEX", ',', 
				eV, 1, fLMOMPT);
		SetEnergyVaryingProperty("./options/mpt/pbwo4-abs.csv", "ABSLENGTH", ',', 
				eV, mm, fLMOMPT);
		BOM = 0;

		std::vector<G4double> ener = {2.6*eV};
		std::vector<G4double> scin = {16000.};

		fLMOMPT->AddProperty("SCINTILLATIONCOMPONENT1",ener,scin);
		fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.005*ns);

		fLMOMPT->AddConstProperty("SCINTILLATIONYIELD",16000./MeV);
		
		fLMOMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
		//fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.*ns);

	} else {

		G4Element* elLi = new G4Element
			(name="Lithium", symbol="Li", iz=3., a=6.94*g/mole);
		G4Element* elMo = new G4Element
			(name="Molybdenum",symbol="Mo", iz=42, a=95.96*g/mole);
		G4Element* elO = new G4Element
			(name="Oxygen", symbol="O", iz=8., a=16.00*g/mole);

		// 280 g, 45mm cube
		G4Material* Li2MoO4 = new G4Material 
			(name="Li2MoO4", density=3.073*g/cm3, 3);

		Li2MoO4->AddElement(elLi, natoms=2);
		Li2MoO4->AddElement(elMo, natoms=1);
		Li2MoO4->AddElement(elO, natoms=4);

		LMOMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Li2MoO4");

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/lmo-rindex.csv", "RINDEX", ',', 
				eV, 1, fLMOMPT);
		SetEnergyVaryingProperty("./options/mpt/lmo-abs.csv", "ABSLENGTH", ',', 
				eV, mm, fLMOMPT);
		BOM = 0;	

		std::vector<G4double> ener = {2.6*eV};//3.*eV};,100.*MeV,10000.*GeV};
		//std::vector<G4double> scin = {1000.};//,100000.,100000.,100000.};
		//std::vector<G4double> scin = {10000.};//,100000.,100000.,100000.};
		std::vector<G4double> scin = {4100.};//,100000.,100000.,100000.};

		fLMOMPT->AddProperty("SCINTILLATIONCOMPONENT1",ener,scin);
		fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.005*ns);

		fLMOMPT->AddConstProperty("SCINTILLATIONYIELD",4100./MeV);
		
		fLMOMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
		//fLMOMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",0.*ns);

		//fLMOMPT->AddProperty("ALPHASCINTILLATIONYIELD",en,ayield); 
		//fLMOMPT->AddProperty("ELECTRONSCINTILLATIONYIELD",en,eyield);
 	
	}

	LMOMaterial->SetMaterialPropertiesTable(fLMOMPT);

	if (! LMOMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLMOMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}

	G4cout << "PRINTING CRYSTAL TABLE\n" << G4endl;
	fLMOMPT->DumpTable();
}



/////////////////
// LD material
/////////////////

void DetectorConstruction::DefineLDM(char mat){
	G4double a,z,density;
	G4String name= "";
	G4String symbol = "";

	G4MaterialPropertiesTable* fLDMPT = new G4MaterialPropertiesTable();

	//G4double trans[4] = {0.1,0.1,0.1,0.1};
	//G4double refl[4] = {0.9,0.9,0.9,0.9}; 

	if(mat == 's'){
		new G4Material(name="Silicon", 
				z=14.0, a=28.0855*g/mole, density=2.33*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Silicon");

		BOM = 0;
		SetEnergyVaryingProperty("./options/mpt/si-rindex.csv", "RINDEX", ',', 
				eV, 1, fLDMPT);
		SetEnergyVaryingProperty("./options/mpt/si-abs.csv", "ABSLENGTH", ',', 
				eV, cm, fLDMPT);
		BOM = 0;	

	} else {
		new G4Material(name="Germanium", 
				z=32.0, a=72.61*g/mole, density=5.323*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Germanium");

		BOM = 0;
		SetEnergyVaryingProperty("./options/mpt/ge-rindex.csv", "RINDEX", ',', 
				eV, 1, fLDMPT);
	
		SetEnergyVaryingProperty("./options/mpt/ge-abs.csv", "ABSLENGTH", ',', 
				eV, cm, fLDMPT);
		BOM = 0;	
	}

	//LDk->AddProperty("REFLECTIVITY", 
	//		energies, refl, numentries);

	LightDetectorMaterial->SetMaterialPropertiesTable(fLDMPT);

	if (! LightDetectorMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLDMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}

	G4cout << "PRINTING LD TABLE\n" << G4endl;
	fLDMPT->DumpTable();
}



/////////////////
// EM material (LMO runs)
/////////////////

void DetectorConstruction::DefineEMM(char mat){
	G4String name= "";
	G4String symbol = "";

	//G4double trans[4] = {0.1,0.1,0.1,0.1};
	//G4double refl[4] = {0.9,0.9,0.9,0.9}; 

	G4MaterialPropertiesTable* EMMPT = new G4MaterialPropertiesTable();

	if(mat == 'c'){

		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Cu");

		BOM = 0;
		SetEnergyVaryingProperty("./options/mpt/cu-rindex.csv", "RINDEX", ',', 
				eV, 1, EMMPT);
		SetEnergyVaryingProperty("./options/mpt/cu-abs.csv", "ABSLENGTH", ',', 
				eV, cm, EMMPT);
		BOM = 0;	

	} else if(mat == 'm') {

		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Al");

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/mysterious-rindex.csv", "RINDEX", ',', 
				eV, 1, EMMPT);
		SetEnergyVaryingProperty("./options/mpt/mysterious-abs.csv", "ABSLENGTH", ',', 
				eV, cm, EMMPT);
		BOM = 0; 

	} else {

		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Al");

		BOM = 0;
		SetEnergyVaryingProperty("./options/mpt/al-rindex.csv", "RINDEX", ',', 
				eV, 1, EMMPT);
		SetEnergyVaryingProperty("./options/mpt/al-abs.csv", "ABSLENGTH", ',', 
				eV, cm, EMMPT);
		BOM = 0; 
	}

	//EMk->AddProperty("TRANSMITTANCE",energies,trans,numentries);
	//EMk->AddProperty("REFLECTIVITY",energies,refl,numentries);

	EMMaterial->SetMaterialPropertiesTable(EMMPT);

	if (!EMMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nEMMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}	

	G4cout << "PRINTING EM TABLE\n" << G4endl;
	EMMPT->DumpTable();
}



/////////////////
// ARC material (LMO runs)
/////////////////

void DetectorConstruction::DefineARCM(char mat){
	G4double density;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	G4MaterialPropertiesTable* ARCMPT = new G4MaterialPropertiesTable();	

	if(mat == '2'){

		G4Element* elSi = G4NistManager::Instance()
			->FindOrBuildElement(14, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		G4Material* SiO2 = new G4Material 
			(name="SiO2", density=2.196*g/cm3, 2);

		SiO2->AddElement(elSi, natoms=1);
		SiO2->AddElement(elO, natoms=2);

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/sio2-rindex.csv", "RINDEX", ',', 
				eV, 1, ARCMPT);
		SetEnergyVaryingProperty("./options/mpt/sio2-abs.csv", "ABSLENGTH", ',', 
				eV, cm, ARCMPT);
		BOM = 0;	

	} else if(mat == 'n'){

		G4Element* elNb = G4NistManager::Instance()
			->FindOrBuildElement(41, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		G4Material* Nb2O5 = new G4Material 
			(name="Nb2O5", density=4.6*g/cm3, 2);

		Nb2O5->AddElement(elNb, natoms=2);
		Nb2O5->AddElement(elO, natoms=5);

		ARCMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Nb2O5");

		BOM = 99;
		SetEnergyVaryingProperty("./options/mpt/nb2o5-rindex.csv", "RINDEX", ',', 
				eV, 1, ARCMPT);
		SetEnergyVaryingProperty("./options/mpt/nb2o5-abs.csv", "ABSLENGTH", ',', 
				eV, mm, ARCMPT);
		BOM = 0;

	}else{

		G4Element* elSi = G4NistManager::Instance()
			->FindOrBuildElement(14, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		G4Material* SiO = new G4Material 
			(name="SiO", density=2.196*g/cm3, 2);

		SiO->AddElement(elSi, natoms=1);
		SiO->AddElement(elO, natoms=1);

		ARCMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("SiO");

		SetEnergyVaryingProperty("./options/mpt/sio-rindex.csv", "RINDEX", ',', 
				eV, 1, ARCMPT);

		//SetEnergyVaryingProperty("./options/mpt/sio-abs.csv", "ABSLENGTH", ',', 
		//		eV, m, ARCMPT);

		//SetEnergyVaryingProperty("./options/mpt/sio-abs.csv", "ABSLENGTH", ',', 
		//		eV, cm, ARCMPT);

		double en[2] = {0.0000000001,100000000000};
		double abs[2] = {0.004675956233*cm,0.004675956233*cm};
		ARCMPT->AddProperty("ABSLENGTH", en, abs, 2);


		/*
		G4int numentries = 4;
		G4double energies[4] = {1.0*eV, 3.0*eV, 4.0*eV,100*GeV};
		ARCMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Vacuum");
		G4double rindex[4] = {1,1,1,1}; 
		G4double abslength[4]={(200)*mm,(200)*mm,(200)*mm,(200)*mm};

		ARCMPT->AddProperty("RINDEX", energies, rindex, numentries);
		ARCMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
		*/

	}

	ARCMaterial->SetMaterialPropertiesTable(ARCMPT);

	G4cout << "PRINTING ARC TABLE\n" << G4endl;
	ARCMPT->DumpTable();
}



/////////////////
// Misc material definitions for CoCalib
/////////////////

void DetectorConstruction::DefineMisc(G4Material** CuMaterial, G4Material** AlMaterial){


	*CuMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4MaterialPropertiesTable* CuMPT = new G4MaterialPropertiesTable();	

	BOM = 0;
	SetEnergyVaryingProperty("./options/mpt/cu-rindex.csv", "RINDEX", ',', 
			eV, 1, CuMPT);
	//SetEnergyVaryingProperty("./options/mpt/cu-abs.csv", "ABSLENGTH", ',', 
	//		eV, cm, CuMPT);
	BOM = 0;	

	double en[2] = {0.0000000001,100000000000};
	double absCu[2] = {0.00000167714884696*cm,0.00000167714884696*cm};
	CuMPT->AddProperty("ABSLENGTH", en, absCu, 2);

	(*CuMaterial)->SetMaterialPropertiesTable(CuMPT);


	*AlMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4MaterialPropertiesTable* AlMPT = new G4MaterialPropertiesTable();	

	double absAl[2] = {0.000000675401864109*cm,0.000000675401864109*cm};
	double rAl[2] = {1120,1120};
	AlMPT->AddProperty("ABSLENGTH", en, absAl, 2);
	AlMPT->AddProperty("RINDEX", en, rAl, 2);

	(*AlMaterial)->SetMaterialPropertiesTable(AlMPT);


	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING Cu TABLE\n" << G4endl;
	CuMPT->DumpTable();
	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING Al TABLE\n" << G4endl;
	AlMPT->DumpTable();
}



/////////////////
// Misc material definitions for the TES module
/////////////////

void DetectorConstruction::DefineTESMisc(G4Material** FrameMaterial, G4Material** LDMaterial, G4Material** ARCMaterial, G4Material** ReflectorMaterial){


	*FrameMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4MaterialPropertiesTable* CuMPT = new G4MaterialPropertiesTable();	

	BOM = 0;
	SetEnergyVaryingProperty("./options/mpt/cu-rindex.csv", "RINDEX", ',', 
			eV, 1, CuMPT);
	//SetEnergyVaryingProperty("./options/mpt/cu-abs.csv", "ABSLENGTH", ',', 
	//		eV, cm, CuMPT);
	BOM = 0;	

	double en[2] = {0.0000000001,100000000000};
	double absCu[2] = {0.00000167714884696*cm,0.00000167714884696*cm};
	CuMPT->AddProperty("ABSLENGTH", en, absCu, 2);

	(*FrameMaterial)->SetMaterialPropertiesTable(CuMPT);
	G4cout << "PROBLEM?" << G4endl;

	// Si LDs
	new G4Material("Silicon", 
			14.0, 28.0855*g/mole, 2.33*g/cm3);

	*LDMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Silicon");
	G4MaterialPropertiesTable* SiMPT = new G4MaterialPropertiesTable();	

	BOM = 0;
	SetEnergyVaryingProperty("./options/mpt/si-rindex.csv", "RINDEX", ',', 
			eV, 1, SiMPT);
//	SetEnergyVaryingProperty("./options/mpt/si-abs.csv", "ABSLENGTH", ',', 
//			eV, cm, SiMPT);
	BOM = 0;	

	double absSi[2] = {0.00131558175025*cm,0.00131558175025*cm};
	SiMPT->AddProperty("ABSLENGTH", en, absSi, 2);

	(*LDMaterial)->SetMaterialPropertiesTable(SiMPT);
	G4cout << "PROLEM?" << G4endl;

	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING LD TABLE\n" << G4endl;
	SiMPT->DumpTable();

	// Si3N4 ARC
	G4Element* elSi = G4NistManager::Instance()
		->FindOrBuildElement(14, false);
	G4Element* elN = G4NistManager::Instance()
		->FindOrBuildElement(7, false);

	G4Material* Si3N4 = new G4Material("Si3N4", 3.17*g/cm3, 2);

	Si3N4->AddElement(elSi, 3);
	Si3N4->AddElement(elN, 4);

	*ARCMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Si3N4");
	G4MaterialPropertiesTable* Si3N4MPT = new G4MaterialPropertiesTable();	

	//double en1[5] = {0.0000000001*eV,0.356*eV,0.5*eV,1*eV,100000000000*eV};
	//double rSi3N4[5] = {1.798,1.927,1.969,2.005,4};
	//Si3N4MPT->AddProperty("RINDEX", en, rSi3N4, 5);

	double en1[8] = {0.0000000001*eV,0.356*eV,0.5*eV,0.734*eV,1*eV,2*eV,5*eV,100000000000*eV};
	double rSi3N4[8] = {1.798,1.927,1.969,1.99,2.005,2.02,2.33,4};
	Si3N4MPT->AddProperty("RINDEX", en1, rSi3N4, 8);

	// extinction coeff < 6*10^-4, 4pik/lamba = a
	// (75.398 * 10^-4) / (400 * 10^7) = 75.398 / (400 * 10^11) = 1.885 * 10^-12 
	//double absSi3N4[2] = {11111*cm,11111*cm};
	double absSi3N4[2] = {0.4080966*cm,0.4080966*cm};
	//double absSi3N4[2] = {0.000000000001885*cm,0.000000000001885*cm};
	Si3N4MPT->AddProperty("ABSLENGTH", en, absSi3N4, 2);

	(*ARCMaterial)->SetMaterialPropertiesTable(Si3N4MPT);

	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING ARC TABLE\n" << G4endl;
	Si3N4MPT->DumpTable();


	/*
	*/
	
	G4cout << "PRONLEM?" << G4endl;



// MISC
	*ReflectorMaterial = *FrameMaterial; 

	G4cout << "PROEM?" << G4endl;

/*
	*AlMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4MaterialPropertiesTable* AlMPT = new G4MaterialPropertiesTable();	

	double absAl[2] = {0.000000675401864109*cm,0.000000675401864109*cm};
double rAl[2] = {1120,1120};
	AlMPT->AddProperty("ABSLENGTH", en, absAl, 2);
	AlMPT->AddProperty("RINDEX", en, rAl, 2);

	(*AlMaterial)->SetMaterialPropertiesTable(AlMPT);

	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING Cu TABLE\n" << G4endl;
	CuMPT->DumpTable();

	G4cout << "\n\n\n" << G4endl;
	G4cout << "PRINTING Al TABLE\n" << G4endl;
	AlMPT->DumpTable();
*/
}



/////////////////
// See above for descriptions 
/////////////////

float* readCSV(const char* filename, char delimiter,int*cols,int*rows){
	FILE *p = fopen(filename,"r");
	char *split = NULL;

	if(p == NULL){
		G4cout << "File unreachable" << G4endl;
		return NULL;
	}

	char te;
	(*rows) = 0;
	(*cols) = 0;

	while((te = fgetc(p)) != EOF){
		//G4cout << "Char: " << te << "\n" << G4endl;

		if((*rows) == 0){
			if(te == delimiter){
				//G4cout << "DELIMITER\n" << G4endl; 
				(*cols)++; 
			}
		}

		if(te == '\n'){
			//G4cout << "NEWLINE\n" << G4endl;
			(*rows)++;
		}
	}

	rewind(p);

	float *arrpnts = (float*)malloc((*rows)*(*cols+1)*sizeof(float*));

	char out[200]; 
	int crows = 0;
	int ccols = 0;

	fgets(out,200,p);


	char * out_s;

	if(BOM == 0){
		out_s = out + 3; // strip BOM
	} else if(BOM == 2){
		out_s = out + 4;
	} else {
		out_s = out;
	}


	do {

		out_s[strlen(out_s)-1] = '\0'; // really can't have line > 200

		split = strtok(out_s,&delimiter);

		while(split){

			if(ccols > (*cols)){
				break;	
			}

			*(arrpnts + ((((*cols)+1)*crows)+ccols)) = atof(split);
			ccols++;
			split = strtok(NULL,&delimiter);
		}

		ccols = 0;
		crows++;

	} while(fgets(out_s,200,p));
	fclose(p);

	G4cout << filename << G4endl;
	for(int i = 0; i < ((*cols)+1)*(*rows); i++){
	   G4cout << "Vals " << i << ": " << *(arrpnts+i) << G4endl;
	 }
	G4cout << "\n\n\n" << G4endl;
	

	return arrpnts;
}



void SetEnergyVaryingProperty(const char* filename, const char* property, const char delimiter, double energyunit, double propertyunit, G4MaterialPropertiesTable* MPT){
	
	int cols, rows, eit, nit;
	float *vals;
	
	vals = readCSV(filename, delimiter, &cols, &rows);

	eit = 0; G4double emen[rows] = {0};
	nit = 0; G4double emrn[rows] = {0};

	for(int i = 0; i < rows*(cols+1); i++){


		//G4cout << "PROPERTY DATA: " << *(vals+i) << G4endl;

		if(i % 2 == 0){// col 0
			emen[eit] = *(vals+i) * energyunit; 
			eit++;
		}

		if((i+1) % 2 == 0){ // col 1
			emrn[nit] = *(vals+i) * propertyunit;
			nit++;
		}
	}


	MPT->AddProperty(property,emen,emrn,rows);
}



void SetEnergyVaryingProperties(const char* filename, const char* property1, const char* property2, const char delimiter, double energyunit, double propertyunit1, double propertyunit2, G4MaterialPropertiesTable* MPT){
	
	int cols, rows, eit, nit, kit;
	float *vals;
	
	vals = readCSV(filename, delimiter, &cols, &rows);

	eit = 0; G4double emen[rows] = {0};
	nit = 0; G4double emrn[rows] = {0};
	kit = 0; G4double emkn[rows] = {0};


	for(int i = 0; i < rows*(cols+1); i++){

		if(i % 3 == 0){// col 0
			emen[eit] = *(vals+i) * energyunit; 
			eit++;
		}

		if((i+1) % 3 == 0){ // col 2
			emkn[kit] = *(vals+i) * propertyunit1;
			kit++;
		}

		if((i+2) % 3 == 0){ // col 1
			emrn[nit] = *(vals+i) * propertyunit2;
			nit++;
		}
	}

	MPT->AddProperty(property1, 
			emen, emrn, rows);
	MPT->AddProperty(property2, 
			emen, emkn, rows);
	
}
