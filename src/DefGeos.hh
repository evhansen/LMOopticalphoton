
static char readExtFile(const char* filename);

// YES, YES to be messed with laterrrrrrr 

G4VPhysicalVolume* DetectorConstruction::DefineCrys(G4double LMOHalfSizes[]){
	auto LMO_box
	= new G4Box("LMO",          // its name
	LMOHalfSizes[0],LMOHalfSizes[1],LMOHalfSizes[2]);    //its size

	G4Material*Mpnt;
	
	char lmos = readExtFile("./options/lmos");

	if(lmos == 'e' || lmos == 's'){ // LMO exists 
		Mpnt = LMOMaterial;
	} else { // LMO doesn't exist
		Mpnt = WorldMaterial;
	}
	
	G4LogicalVolume* fLMO_LV_primary
	= new G4LogicalVolume(
	LMO_box,         // its solid
	Mpnt,    // its material
	"LMO_primary");          // its name

	if(lmos == 'e' || lmos == 's'){ // LMO exists 
		G4VisAttributes* LMOVA = new G4VisAttributes();
		LMOVA->SetColor(1.0,1.0,0.0,0.3);
		LMOVA->SetForceSolid(true);
		fLMO_LV_primary->SetVisAttributes(LMOVA);
	}

	G4VPhysicalVolume* fLMO1 = new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(0,0,0),  // place off center eventually
	fLMO_LV_primary,          // its logical volume
	"LMO1",          // its name
	worldLV,                // its mother  volume
	false,            // no boolean operation
	1,                // copy number
	fCheckOverlaps);  // checking overlaps
	
	return fLMO1;
}

void DetectorConstruction::DefineIV(G4double IVHalfSizes[],G4double IVPos[],
G4VPhysicalVolume** IV1,G4VPhysicalVolume** IV2,G4VPhysicalVolume** IV3,G4VPhysicalVolume** IV4,G4VPhysicalVolume** IV5,G4VPhysicalVolume** IV6){
	
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");

		
	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);

	// for IV between LMO and LD

	G4VSolid* XInterimVacuum = new G4Box("LDInterimVacuum",
		IVHalfSizes[0],IVHalfSizes[1], (IVHalfSizes[2]>0) ? IVHalfSizes[2] : 1*mm ); 

	// yes, I'm having a hernia
		


	G4LogicalVolume* LDIVLV
	= new G4LogicalVolume(
	XInterimVacuum,     // its solid
	WorldMaterial,  // its material
	"XInterimVacuumLV");   // its name

	G4VisAttributes* IVVA = new G4VisAttributes();
	IVVA->SetColor(1.0,0.0,1.0,0.15);
	IVVA->SetForceSolid(true);
	
	LDIVLV->SetVisAttributes(IVVA);
	
	// for IV between LMO and EM
	G4VSolid* EMInterimVacuum = new G4Box("EMInterimVacuum",
		IVHalfSizes[0],IVHalfSizes[1], (IVHalfSizes[3]>0) ? IVHalfSizes[3] : 1*mm );

	
	G4LogicalVolume* EMIVLV
	= new G4LogicalVolume(
	EMInterimVacuum,     // its solid
	WorldMaterial,  // its material
	"EMInterimVacuumLV");   // its name

	EMIVLV->SetVisAttributes(IVVA);

	// I know that this is ugly .... will change after it functions as intended
	
	if (face1 == 'e'){
			
		if(IVHalfSizes[3] > 0){
			*IV1 = new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(IVPos[0],
				IVPos[1],IVPos[2]),
			EMIVLV,          // its logical volume
			"IV1T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			30,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV1 = new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(IVPos[0],
				IVPos[1],IVPos[2]),
			LDIVLV,          // its logical volume
			"IV1T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			30,                // copy number
			fCheckOverlaps);  // checking overlaps					
		}	
	}

	if (face2 == 'e'){
		
		if(IVHalfSizes[3] > 0){
			*IV2 = new G4PVPlacement(
			0,                // no rotation 
			G4ThreeVector(IVPos[3],
				IVPos[4],IVPos[5]),
			EMIVLV,          // its logical volume
			"IV1B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			31,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV2 = new G4PVPlacement(
			0,                // no rotation 
			G4ThreeVector(IVPos[3],
				IVPos[4],IVPos[5]),
			LDIVLV,          // its logical volume
			"IV1B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			31,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	}

	if (face3 == 'e'){
		
		if(IVHalfSizes[3] > 0){
			*IV3 = new G4PVPlacement(
			t1,                // 
			G4ThreeVector(IVPos[6],
				IVPos[7],IVPos[8]),
			EMIVLV,          // its logical volume
			"IV2B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			32,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV3 = new G4PVPlacement(
			t1,                // 
			G4ThreeVector(IVPos[6],
				IVPos[7],IVPos[8]),
			LDIVLV,          // its logical volume
			"IV2B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			32,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	}

	if (face4 == 'e'){
		
		if(IVHalfSizes[3] > 0){
			*IV4 = new G4PVPlacement(
			t1,                // 
			G4ThreeVector(IVPos[9],
				IVPos[10],IVPos[11]),
			EMIVLV,          // its logical volume
			"IV3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			33,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV4 = new G4PVPlacement(
			t1,                // 
			G4ThreeVector(IVPos[9],
				IVPos[10],IVPos[11]),
			LDIVLV,          // its logical volume
			"IV3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			33,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	}

	if (face5 == 'e'){
		
		if(IVHalfSizes[3] > 0){
			*IV5 = new G4PVPlacement(
			p1,                // 
			G4ThreeVector(IVPos[12],
				IVPos[13],IVPos[14]),
			EMIVLV,          // its logical volume
			"IV3B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			34,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV5 = new G4PVPlacement(
			p1,                // 
			G4ThreeVector(IVPos[12],
				IVPos[13],IVPos[14]),
			LDIVLV,          // its logical volume
			"IV3B",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			34,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	}

	if (face6 == 'e'){
		
		if(IVHalfSizes[3] > 0){
			*IV6 = new G4PVPlacement(
			p1,                // 
			G4ThreeVector(IVPos[15],
				IVPos[16],IVPos[17]),
			EMIVLV,          // its logical volume
			"IV3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			35,                // copy number
			fCheckOverlaps);  // checking overlaps
		}
	} else {
		
		if(IVHalfSizes[2] > 0){
			*IV6 = new G4PVPlacement(
			p1,                // 
			G4ThreeVector(IVPos[15],
				IVPos[16],IVPos[17]),
			LDIVLV,          // its logical volume
			"IV3T",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			35,                // copy number
			fCheckOverlaps);  // checking overlaps
		}	
	}
}


void DetectorConstruction::DefineLD(G4double LDHalfSizes[],G4double LDPos[],
G4VPhysicalVolume** LD1,G4VPhysicalVolume** LD2,G4VPhysicalVolume** LD3,G4VPhysicalVolume** LD4,G4VPhysicalVolume** LD5,G4VPhysicalVolume** LD6){
	
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");

		
	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);
	
	auto LightDetector
	= new G4Box("LightDetector",     // its name
	LDHalfSizes[0],LDHalfSizes[1],LDHalfSizes[2]); // its size

	G4LogicalVolume* LightDetectorLV
	= new G4LogicalVolume(
	LightDetector,     // its solid
	LightDetectorMaterial,  // its material
	"LDLV");   // its name

	auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	simpleBoxVisAtt->SetVisibility(true);
	LightDetectorLV->SetVisAttributes(simpleBoxVisAtt);

	if(face1 == 'l'){
		*LD1 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(LDPos[0],
			LDPos[1],LDPos[2]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD1T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		10,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	if(face2 == 'l'){
		*LD2 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(LDPos[3],
			LDPos[4],LDPos[5]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD1B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		11,                // copy number
		fCheckOverlaps);  // checking overlaps
	}

	if(face3 == 'l'){
		*LD3 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(LDPos[6],
			LDPos[7],LDPos[8]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD2T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		12,                // copy number
		fCheckOverlaps);  // checking overlaps
	}

	if(face4 == 'l'){
		*LD4 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(LDPos[9],
			LDPos[10],LDPos[11]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD2B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		13,                // copy number
		fCheckOverlaps);  // checking overlaps
	}

	if(face5 == 'l'){
		*LD5 = new G4PVPlacement(
		p1,                // no rotation
		G4ThreeVector(LDPos[12],
			LDPos[13],LDPos[14]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD3T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		14,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	if(face6 == 'l'){
		*LD6 = new G4PVPlacement(
		p1,                // no rotation
		G4ThreeVector(LDPos[15],
			LDPos[16],LDPos[17]),  // at (0,0,0)
		LightDetectorLV,          // its logical volume
		"LD3B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		15,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
}

void DetectorConstruction::DefineEM(G4double EMHalfSizes[],G4double EMPos[],
G4VPhysicalVolume** EM1,G4VPhysicalVolume** EM2,G4VPhysicalVolume** EM3,G4VPhysicalVolume** EM4,G4VPhysicalVolume** EM5,G4VPhysicalVolume** EM6){
	
	char face1 = readExtFile("./options/face1");
	char face2 = readExtFile("./options/face2");
	char face3 = readExtFile("./options/face3");
	char face4 = readExtFile("./options/face4");
	char face5 = readExtFile("./options/face5");
	char face6 = readExtFile("./options/face6");

	
	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);
	
	auto EM
	= new G4Box("EncapsulatingMaterial",
	EMHalfSizes[0],EMHalfSizes[1],EMHalfSizes[2]);

	G4LogicalVolume* EMLV = new G4LogicalVolume(
	EM,     // its solid
	EMMaterial,  // its material
	"EMLV");   // its name

	G4VisAttributes* EMVA = new G4VisAttributes();
	EMVA->SetColor(0.0,1.0,0.0,1.0);
	EMVA->SetForceLineSegmentsPerCircle(4);
	EMLV->SetVisAttributes(EMVA);

	if(face1 == 'e'){
		*EM1 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(EMPos[0],
			EMPos[1],EMPos[2]),
		EMLV,          
		"EM1T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		20,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	if(face2 == 'e'){
		*EM2 = new G4PVPlacement(
		0,                // no rotation
		G4ThreeVector(EMPos[3],
			EMPos[4],EMPos[5]),  
		EMLV,          // its logical volume
		"EM1B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		21,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	if(face3 == 'e'){
		*EM3 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(EMPos[6],
			EMPos[7],EMPos[8]),
		EMLV,          // its logical volume
		"EM2T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		22,                // copy number
		fCheckOverlaps);  // checking overlaps
	}

	if(face4 == 'e'){
		*EM4 = new G4PVPlacement(
		t1,                //
		G4ThreeVector(EMPos[9],
			EMPos[10],EMPos[11]), 
		EMLV,          // its logical volume
		"EM2B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		23,                // copy number
		fCheckOverlaps);  // checking overlaps
	}

	if(face5 == 'e'){
		*EM5 = new G4PVPlacement(
		p1,                //
		G4ThreeVector(EMPos[12],
			EMPos[13],EMPos[14]), 
		EMLV,          // its logical volume
		"EM3T",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		24,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
	
	if(face6 == 'e'){
		*EM6 = new G4PVPlacement(
		p1,                //
		G4ThreeVector(EMPos[15],
			EMPos[16],EMPos[17]), 
		EMLV,          // its logical volume
		"EM3B",    // its name
		worldLV,          // its mother  volume
		false,            // no boolean operation
		25,                // copy number
		fCheckOverlaps);  // checking overlaps
	}
}

void DetectorConstruction::DefineARC(G4double ARCHalfSizes[],G4double ARCPos[],
G4VPhysicalVolume** ARC1,G4VPhysicalVolume** ARC2,G4VPhysicalVolume** ARC3,G4VPhysicalVolume** ARC4,G4VPhysicalVolume** ARC5,G4VPhysicalVolume** ARC6){
	
		
	G4RotationMatrix* t1 = new G4RotationMatrix();
	t1->rotateY(90.*deg);
	t1->rotateX(90.*deg);

	G4RotationMatrix* p1 = new G4RotationMatrix();
	p1->rotateY(90.*deg);
	p1->rotateZ(90.*deg);
	
	auto ARC = new G4Box("Antirefcoat",
	ARCHalfSizes[0],ARCHalfSizes[1],ARCHalfSizes[2]);

	G4LogicalVolume* ARCLV
	= new G4LogicalVolume(
	ARC,     // its solid
	ARCMaterial,  // its material
	"ARCLV");   // its name

	G4VisAttributes* ARCVA = new G4VisAttributes();
	ARCVA->SetColor(1.0,1.0,1.0,0.5);
	ARCVA->SetForceSolid(true);
	ARCLV->SetVisAttributes(ARCVA);

	*ARC1 = new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(ARCPos[0],
		ARCPos[1],ARCPos[2]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC1T",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	40,                // copy number
	fCheckOverlaps);  // checking overlaps

	*ARC2 = new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(ARCPos[3],
		ARCPos[4],ARCPos[5]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC1B",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	41,                // copy number
	fCheckOverlaps);  // checking overlaps

	*ARC3 = new G4PVPlacement(
	t1,                //
	G4ThreeVector(ARCPos[6],
		ARCPos[7],ARCPos[8]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC2T",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	42,                // copy number
	fCheckOverlaps);  // checking overlaps

	*ARC4 = new G4PVPlacement(
	t1,                //
	G4ThreeVector(ARCPos[9],
		ARCPos[10],ARCPos[11]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC2B",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	43,                // copy number
	fCheckOverlaps);  // checking overlaps

	*ARC5 = new G4PVPlacement(
	p1,                // no rotation
	G4ThreeVector(ARCPos[12],
		ARCPos[13],ARCPos[14]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC3T",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	44,                // copy number
	fCheckOverlaps);  // checking overlaps
				
	*ARC6 = new G4PVPlacement(
	p1,                // no rotation
	G4ThreeVector(ARCPos[15],
		ARCPos[16],ARCPos[17]),  // at (0,0,0)
	ARCLV,          // its logical volume
	"ARC3B",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	45,                // copy number
	fCheckOverlaps);  // checking overlaps
}



/*void DetectorConstruction::DefineSurfaces(G4VPhysicalVolume* EM1,G4VPhysicalVolume* EM2,G4VPhysicalVolume* EM3,G4VPhysicalVolume* EM4,G4VPhysicalVolume* EM5,G4VPhysicalVolume* EM6,G4VPhysicalVolume* ARC1,G4VPhysicalVolume* ARC2,G4VPhysicalVolume* ARC3,G4VPhysicalVolume* ARC4,G4VPhysicalVolume* ARC5,G4VPhysicalVolume* ARC6){
	
}*/



