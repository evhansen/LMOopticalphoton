#include "LightDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4PrimaryParticle.hh"
#include "G4RunManager.hh"

#include <stddef.h>


void populateHit(const G4Step* step, G4ParticleDefinition*& pdef,
	G4double&edep, G4double&stepLength, G4double&energy,
	G4ThreeVector&pos, G4ThreeVector&momdir,
	G4int&pdg, G4int&copynum, G4int&eventid, G4int&trackid, G4int&parentid,
	std::string&x, std::string&y, std::string&z, 
	std::string&px, std::string&py, std::string&pz);
void chkCols(std::string line,std::string&x,std::string&y,std::string&z,
	std::string&px,std::string&py,std::string&pz);
std::string getCol(std::string line, unsigned int col);
std::string readExtFileLine(const char* filename);
std::string readLast(const char* fname, G4int eid);



LightDetectorSD::LightDetectorSD(const G4String& name,
	const G4String& hitsCollectionName,G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);

}

LightDetectorSD::~LightDetectorSD()
{
	std::string fname = "./data/" + readExtFileLine("./options/name");
	remove(fname.c_str());
}


void LightDetectorSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection // LightDetectorHitsCollection
	  //= new G4THitsCollection<LightDetectorHit>;
	  = new LightDetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

}


G4bool LightDetectorSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{

	G4double edep, stepLength, energy;
	G4ParticleDefinition* pdef;
	G4ThreeVector pos, momdir;
	G4int pdg, copynum, eventid, trackid, parentid;
	std::string px,py,pz,x,y,z;


	populateHit(step, pdef,
		edep, stepLength, energy,
		pos, momdir,
		pdg, copynum, eventid, trackid, parentid,
		x,y,z,px,py,pz);


	if(edep == 0.){return false;}
	if (stepLength == 0.) {return false;}
 
 
	// populate hit
	LightDetectorHit* OPhit = new LightDetectorHit();
 
	OPhit->SetEdep(edep);
	OPhit->SetEnergy(energy);
  
	OPhit->SetPDef(pdef);
  
	OPhit->SetPos(pos);
	OPhit->SetfMom(momdir);

	OPhit->SetPDG(pdg);
	OPhit->SetPhysVolNum(copynum);
	
	OPhit->SetEventID(eventid);
	OPhit->SetTrackID(trackid);
	OPhit->SetParentID(parentid);

	
	// convert units m -> mm
	OPhit->SetiPos(G4ThreeVector(stof(x)*1000,stof(y)*1000,stof(z)*1000));	

	OPhit->SetiMom(G4ThreeVector(stof(px),stof(py),stof(pz)));

	fHitsCollection->insert(OPhit);


	// this even necessary? No idea
	G4int snum = step->GetNumberOfSecondariesInCurrentStep();
	if(snum == 0){return false;} //want this?
	
	const std::vector<const G4Track*>* secs = step->GetSecondaryInCurrentStep();

	LightDetectorHit** LDH = (LightDetectorHit**)malloc(snum*sizeof(LightDetectorHit*)); 

	//TODO: add primary vertices to event for secondaries
	//G4Event* sevent = G4RunManager::GetRunManager()->GetCurrentEvent();
  
	for(G4int i = 0; i < snum; i++){
	
		const G4Step* sstep = (*secs)[i]->GetStep();
		if(sstep == NULL){continue;}
	
		// track = (*secs)[i];
		populateHit(sstep, pdef,
			edep, stepLength, energy,
			pos, momdir,
			pdg, copynum, eventid, trackid, parentid,
			x,y,z,px,py,pz);


		if(edep == 0.){continue;}  
		if (stepLength == 0.) {continue;}

		LDH[i] = new LightDetectorHit();
 
		LDH[i]->SetEdep(edep);
	  	LDH[i]->SetEnergy(energy);

		LDH[i]->SetPDef(pdef);
	  
		LDH[i]->SetPos(pos);
		LDH[i]->SetfMom(momdir);

	  	LDH[i]->SetPDG(pdg);
		LDH[i]->SetPhysVolNum(copynum);

		LDH[i]->SetEventID(eventid);
		LDH[i]->SetTrackID(trackid);
		LDH[i]->SetParentID(parentid);


		// convert units m -> mm
		LDH[i]->SetiPos(G4ThreeVector(stof(x)*1000,stof(y)*1000,stof(z)*1000));

		LDH[i]->SetiMom(G4ThreeVector(stof(px),stof(py),stof(pz)));

		fHitsCollection->insert(LDH[i]);
	}



  return true;
}


void LightDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl
       << "-------->Hits Collection: in this event they are " << nofHits
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}




void populateHit(const G4Step* step,G4ParticleDefinition*& pdef,
		G4double&edep, G4double&stepLength, G4double&energy,
		G4ThreeVector&pos, G4ThreeVector&momdir,
		G4int&pdg, G4int&copynum, G4int&eventid, G4int&trackid, G4int&parentid,
		std::string&x, std::string&y, std::string&z, 
		std::string&px, std::string&py, std::string&pz){

	if(step == NULL){return;}
	
	edep = step->GetTotalEnergyDeposit();
	stepLength = step->GetStepLength();

	if(edep == 0.){return;}
	if (stepLength == 0.) {return;}
	
	G4Track* track = step->GetTrack();
	
	energy = track->GetDynamicParticle()->GetTotalEnergy();

	pdef = track->GetDefinition();

	G4StepPoint* thePrePoint = step->GetPreStepPoint();
	G4StepPoint* thePostPoint = step->GetPostStepPoint();
	
	pos = (thePrePoint->GetPosition() + thePostPoint->GetPosition()) / 2.;
	momdir = track->GetDynamicParticle()->GetMomentumDirection();

	pdg = track->GetDynamicParticle()->GetPDGcode(); 
	copynum = thePrePoint->GetTouchable()->GetCopyNumber();

	eventid = -1;
	
	if(G4RunManager::GetRunManager()->GetCurrentEvent()){
		eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	}

	trackid = track->GetTrackID();
	parentid = track->GetParentID();


	// checking file
	static std::string fname = "./data/" + readExtFileLine("./options/name");

	std::string line;

	x = "0";
	y = "0";
	z = "0";
	px = "0";
	py = "0";
	pz = "0";


	if((parentid == 0) && (eventid > -1)){
	
		line = readLast(fname.c_str(),eventid);

		if(line != "no"){

			// 0 is eventid
			px = getCol(line,1);
			py = getCol(line,2);
			pz = getCol(line,3);
			x = getCol(line,4);
			y = getCol(line,5);
			z = getCol(line,6);

			chkCols(line,x,y,z,px,py,pz);
		}
	}

}


void chkCols(std::string line,std::string&x,std::string&y,std::string&z,
		std::string&px,std::string&py,std::string&pz){

	if(px == "no" || py == "no" || pz == "no" || x == "no" || y == "no" || z == "no"){ 
	
		G4cout << "******************************* LINE ISSUE **************************\n" 
			<< line << "\n"; 
		G4cout << px <<":"<< py <<":"<< pz<<"\n";
		G4cout << x <<":"<< y <<":"<< z<<"\n";
		
		px = "0";
		py = "0";
		pz = "0";
		x = "0";
		y = "0";
		z = "0";
	}
}


std::string getCol(std::string line, unsigned int col){

	std::string sstr;

	int b = 0; int e = line.find(",",0);
	
	if(col == 0){return line.substr(b,e);}

	unsigned int iter = 1;
	
	while(b < line.length()){

		b = e+1; e = line.find(",",b); 
		
		if(col == iter){return line.substr(b,e-b);}
		
		iter++;

		if( b > e || b == e ) {break;}
	}

	return "no";
}


std::string readLast(const char* fname, G4int eid){
	std::ifstream p;
	std::string line;
	std::string event;

	p.open(fname);

	char ch = '0';
	p.seekg(-1,std::ios_base::end);

	int loc = -1;

	while((p.tellg() > 1)) {
	
		while(ch != '\n'){
			p.get(ch); 
			--loc; 
			p.seekg(loc,std::ios_base::end); 
			p.get(ch);
		}

		getline(p,line);
 
	 	event = line.substr(0,line.find(","));

		if(atoi(event.c_str()) == eid){ return line;}

		--loc;		
		p.seekg(loc,std::ios_base::end);
		p.get(ch);
	}
	

	p.close();

	return "no";
}





