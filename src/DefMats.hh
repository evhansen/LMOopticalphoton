
static int BOM = 0; // 0 if BOM needs to be removed, else do not remove.


static char readExtFile(const char* filename);
float* readCSV(const char* filename, char delimiter,int*cols,int*rows);


//TODO: fix up these

void DetectorConstruction::DefineWorldM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {1.0*eV, 3.0*eV, 4.0*eV,100*GeV};
	G4double a,z,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;


	new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);

	G4MaterialPropertiesTable* fWorldMPT = new G4MaterialPropertiesTable();

	G4double rindex[4] = {1,1,1,1}; 
	G4double abslength[4]={200*mm,200*mm,200*mm,200*mm};

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
}




void DetectorConstruction::DefineCrysM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {1.0*eV, 3.0*eV, 4.0*eV,100*GeV};
	G4double a,z,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	int cols, rows, eit, nit, kit;
	float *vals;

	G4double abslength[4], rindex[4];


	G4MaterialPropertiesTable* fLMOMPT = new G4MaterialPropertiesTable();	


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


		abslength[0] = 5*mm; abslength[1] = 5*mm;
		abslength[2] = 5*mm; abslength[3] = 80*cm;

		if(readExtFile("./options/teo2rndx") == 'o'){

			vals = readCSV("./options/teo2-norindex.csv",',',&cols,&rows);

		} else {
			vals = readCSV("./options/teo2-nerindex.csv",',',&cols,&rows);
		}	

		eit = 0; G4double cuben[rows] = {0};
		nit = 0; G4double cubrn[rows] = {0};
		//kit = 0; G4double cubkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 2 == 0){// col 0
				cuben[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 2 == 0){ // col 1
				cubrn[nit] = *(vals+i);
				nit++;
			}
		}


		fLMOMPT->AddProperty("RINDEX", 
				cuben, cubrn, rows);
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

		abslength[0] = 100*mm; abslength[1] = 100*mm;
		abslength[2] = 100*mm; abslength[3] = 100*mm;
		rindex[0] = 1.44;	rindex[1] = rindex[0];
		rindex[2] = rindex[0];	rindex[3] = rindex[0];

		fLMOMPT->AddProperty("RINDEX", energies, rindex, numentries);
		fLMOMPT->AddProperty("REALRINDEX", energies, rindex, numentries);
		
		std::vector<G4double> en = {1.*eV,100.*GeV };
		std::vector<G4double> eyield = {250.,250.};// 0.65keV/MeV ~2.6
		std::vector<G4double> ayield = {38.,38.};// 0.1keV/Mev ~2.6

		fLMOMPT->AddConstProperty("SCINTILLATIONYIELD",100000./MeV); // gammas

		//fLMOMPT->AddProperty("ALPHASCINTILLATIONYIELD",en,ayield); 
		//fLMOMPT->AddProperty("ELECTRONSCINTILLATIONYIELD",en,eyield);
 	}

	fLMOMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
	LMOMaterial->SetMaterialPropertiesTable(fLMOMPT);


	if (! LMOMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLMOMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}

}




void DetectorConstruction::DefineLDM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {2.0*eV, 3.0*eV, 4.0*eV,100*MeV};
	G4double a,z,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	int cols, rows, eit, nit, kit;
	float *vals;

	G4double abslength[4];



	G4MaterialPropertiesTable* fLDMPT = new G4MaterialPropertiesTable();

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

	if(mat == 's'){
		new G4Material(name="Silicon", 
				z=14.0, a=28.0855*g/mole, density=2.33*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Silicon");


		vals = readCSV("./options/al-rindex.csv",',',&cols,&rows);

		eit = 0; G4double lden[rows] = {0};
		nit = 0; G4double ldrn[rows] = {0};
		kit = 0; G4double ldkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				lden[eit] = *(vals+i) * eV; 
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
		  Rrindex[0] = 3.906;	Rrindex[1] = 5.222;
		  Rrindex[2] = 5.01; 	Rrindex[3] = 1.01;	
		  Irindex[0] = 0.022; 	Irindex[1] = 0.269;
		  Irindex[2] = 3.586; 	Irindex[3] = 2.909;

		abslength[0] = pow(10,2)*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		// TODO: FIX THE FILE RINDEX, REEEEEEEEEEEEEEE

		fLDMPT->AddProperty("RINDEX", 
			energies, Rrindex, numentries);
		fLDMPT->AddProperty("REALRINDEX", 
			energies, Rrindex, numentries);
		fLDMPT->AddProperty("IMAGINARYRINDEX", 	
			energies, Irindex, numentries);
		LDk->AddProperty("IMAGINARYRINDEX", 
			energies,Irindex, numentries);
		
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

	} else {

		new G4Material(name="Germanium", 
				z=32.0, a=72.61*g/mole, density=5.323*g/cm3);

		LightDetectorMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("Germanium");


		vals = readCSV("./options/ge-rindex.csv",',',&cols,&rows);

		eit = 0; G4double lden[rows] = {0};
		nit = 0; G4double ldrn[rows] = {0};
		kit = 0; G4double ldkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				lden[eit] =*(vals+i); 
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

	//LDk->AddProperty("TRANSMITTANCE", 
	//		energies, trans, numentries);
	LDk->AddProperty("REFLECTIVITY", 
			energies, refl, numentries);

	LightDetectorMaterial->SetMaterialPropertiesTable(fLDMPT);

	if (! LightDetectorMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nLDMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}
}




void DetectorConstruction::DefineEMM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {2.0*eV, 3.0*eV, 4.0*eV,100*MeV};
	G4double a,z,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	int cols, rows, eit, nit, kit;
	float *vals;

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


	G4double abslength[4];


	G4MaterialPropertiesTable* EMMPT = new G4MaterialPropertiesTable();

	if(mat == 'c'){

		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Cu");

		vals = readCSV("./options/cu-rindex.csv",',',&cols,&rows);

		eit = 0; G4double emen[rows] = {0};
		nit = 0; G4double emrn[rows] = {0};
		kit = 0; G4double emkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				emen[eit] = *(vals+i) * eV; 
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

		abslength[0] = 10*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];


		EMMPT->AddProperty("REALRINDEX",emen,emrn,rows);
		EMMPT->AddProperty("RINDEX",emen,emrn,rows);
		EMMPT->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
		EMMPT->AddProperty("ABSLENGTH",energies, abslength, numentries);
		EMk->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);

	} else if(mat == 'm') {

		BOM = 99;
		
		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Al");
	
		vals = readCSV("./options/mysterious-rindex.csv",',',&cols,&rows);	

		eit = 0; G4double emen[rows] = {0};
		nit = 0; G4double emrn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 2 == 0){// col 0
				emen[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 2 == 0){ // col 1
				emrn[nit] = *(vals+i);
				nit++;
			}
		}

		EMMPT->AddProperty("REALRINDEX",emen,emrn,rows);
		EMMPT->AddProperty("RINDEX",emen,emrn,rows);

		cols = 0; rows = 0;	

		vals = readCSV("./options/mysterious-abs.csv",',',&cols,&rows);

		eit = 0; G4double emabsen[rows] = {0};
		nit = 0; G4double emabs[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 2 == 0){// col 0
				emabsen[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 2 == 0){ // col 1
				emabs[nit] = *(vals+i) * mm;
				nit++;
			}
		}

		EMMPT->AddProperty("ABSLENGTH",emabsen, emabs, rows);

		BOM = 0; // reset to default 

	} else {

		EMMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("G4_Al");


		vals = readCSV("./options/al-rindex.csv",',',&cols,&rows);
		eit = 0; G4double emen[rows] = {0};
		nit = 0; G4double emrn[rows] = {0};
		kit = 0; G4double emkn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				emen[eit] = *(vals+i) * eV; 
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

		abslength[0] = pow(10,3)*pow(10,-6)*m;abslength[1] = abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		EMMPT->AddProperty("REALRINDEX",emen,emrn,rows);
		EMMPT->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
		EMMPT->AddProperty("RINDEX",emen,emrn,rows);
		EMMPT->AddProperty("ABSLENGTH",energies, abslength, numentries);
		EMk->AddProperty("IMAGINARYRINDEX",emen,emkn,rows);
	}


	//EMk->AddProperty("TRANSMITTANCE",energies,trans,numentries);
	EMk->AddProperty("REFLECTIVITY",energies,refl,numentries);


	EMMaterial->SetMaterialPropertiesTable(EMMPT);


	if (!EMMaterial) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.\nEMMAT\n";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}	
}




void DetectorConstruction::DefineARCM(char mat){
	G4int numentries = 4;

	G4double energies[4] = {2.0*eV, 3.0*eV, 4.0*eV,100*MeV};
	G4double a,z,density, iz;
	G4String name= "";
	G4String symbol = "";
	G4int natoms;

	int cols, rows, eit, nit, kit;
	float *vals;

	G4double abslength[4];


	G4MaterialPropertiesTable* ARCMPT= new G4MaterialPropertiesTable();	

	if(true){
		G4Element* elSi = G4NistManager::Instance()
			->FindOrBuildElement(14, false);
		G4Element* elO = G4NistManager::Instance()
			->FindOrBuildElement(8, false);

		G4Material* SiO2 = new G4Material 
			(name="SiO2", density=2.196*g/cm3, 2);

		SiO2->AddElement(elSi, natoms=1);
		SiO2->AddElement(elO, natoms=2);

		ARCMaterial = G4NistManager::Instance()
			->FindOrBuildMaterial("SiO2");

		abslength[0] = 1/(383.94*cm);abslength[1]= abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		vals = readCSV("./options/sio2-rindex.csv",',',&cols,&rows);

		eit = 0; G4double arcen[rows] = {0};
		nit = 0; G4double arcrn[rows] = {0};
		kit = 0; G4double arckn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				arcen[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 3 == 0){ // col 2
				arckn[kit] = *(vals+i);
				kit++;
			}

			if((i+2) % 3 == 0){ // col 1
				arcrn[nit] = *(vals+i);
				nit++;
			}
		}

		ARCMPT->AddProperty("RINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("REALRINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
		ARCk->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
	} else if(false){

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

		abslength[0] = 1/(0.004676*cm);abslength[1]= abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		vals = readCSV("./options/sio-rindex.csv",',',&cols,&rows);

		eit = 0; G4double arcen[rows] = {0};
		nit = 0; G4double arcrn[rows] = {0};
		kit = 0; G4double arckn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				arcen[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 3 == 0){ // col 2
				arckn[kit] = *(vals+i);
				kit++;
			}

			if((i+2) % 3 == 0){ // col 1
				arcrn[nit] = *(vals+i);
				nit++;
			}
		}

		ARCMPT->AddProperty("RINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("REALRINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
		ARCk->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
	}else{

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


		abslength[0] = 1/(4.676*cm);abslength[1]= abslength[0];
		abslength[2] = abslength[0];abslength[3] = abslength[0];

		vals = readCSV("./options/nb2o5-rindex.csv",',',&cols,&rows);

		eit = 0; G4double arcen[rows] = {0};
		nit = 0; G4double arcrn[rows] = {0};
		kit = 0; G4double arckn[rows] = {0};

		for(int i = 0; i < rows*(cols+1); i++){

			if(i % 3 == 0){// col 0
				arcen[eit] = *(vals+i) * eV; 
				eit++;
			}

			if((i+1) % 3 == 0){ // col 2
				arckn[kit] = *(vals+i);
				kit++;
			}

			if((i+2) % 3 == 0){ // col 1
				arcrn[nit] = *(vals+i);
				nit++;
			}
		}

		ARCMPT->AddProperty("RINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("REALRINDEX", 
				arcen, arcrn, rows);
		ARCMPT->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
		ARCk->AddProperty("IMAGINARYRINDEX", 
				arcen, arckn, rows);
	}

	ARCMPT->AddProperty("ABSLENGTH", energies, abslength, numentries);
	ARCMaterial->SetMaterialPropertiesTable(ARCMPT);
}








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
		if((*rows) == 0){
			if(te == delimiter){
				(*cols)++; 
			}
		}

		if(te == '\n'){
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

	/*G4cout << filename << G4endl;
	for(int i = 0; i < ((*cols)+1)*(*rows); i++){
	   G4cout << "Vals " << i << ": " << *(arrpnts+i) << G4endl;
	 }
	G4cout << "\n\n\n" << G4endl;
	*/

	return arrpnts;
}
