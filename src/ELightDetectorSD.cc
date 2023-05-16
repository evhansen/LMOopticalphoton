#include "ELightDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4PrimaryParticle.hh"
#include "G4RunManager.hh"

#include <stddef.h>


/////////////////
// Populate the hits
/////////////////
void populateHit(const G4Step* step, G4ParticleDefinition*& pdef,
	G4double&edep, G4double&stepLength, G4double&energy,
	G4ThreeVector&pos, G4ThreeVector&momdir,
	G4int&pdg, G4int&copynum, G4int&eventid, G4int&trackid, G4int&parentid,
	std::string&x, std::string&y, std::string&z, 
	std::string&px, std::string&py, std::string&pz,
	G4String& particlename, G4String& processname);

std::string readExtFileLine(const char* filename);



ELightDetectorSD::ELightDetectorSD(const G4String& name,
	const G4String& hitsCollectionName,G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);

}

ELightDetectorSD::~ELightDetectorSD()
{
	//std::string fname = "./data/" + readExtFileLine("./options/name");
	//remove(fname.c_str());
}


void ELightDetectorSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
	  = new LightDetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

}


G4bool ELightDetectorSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{

	G4double edep, stepLength, energy;
	G4ParticleDefinition* pdef;
	G4ThreeVector pos, momdir;
	G4int pdg, copynum, eventid, trackid, parentid;
	G4String partname, procname;
	std::string px,py,pz,x,y,z;


	populateHit(step, pdef,
		edep, stepLength, energy,
		pos, momdir,
		pdg, copynum, eventid, trackid, parentid,
		x,y,z,px,py,pz,
		partname,procname);


	//if(edep == 0.){return false;}
	if (stepLength == 0.) {return false;}
 
	/////////////////
	// Create the hit
	/////////////////
	LightDetectorHit* OPhit = new LightDetectorHit();

	edep = energy;
 
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

	partname = step->GetTrack()->GetParticleDefinition()->GetParticleName();
	
	OPhit->SetParticleName(partname);
	OPhit->SetProcessName(procname);

	
	// convert units m -> mm
	//OPhit->SetiPos(G4ThreeVector(stof(x)*1000,stof(y)*1000,stof(z)*1000));	

	//OPhit->SetiMom(G4ThreeVector(stof(px),stof(py),stof(pz)));

	fHitsCollection->insert(OPhit);

	step->GetTrack()->SetTrackStatus(fStopAndKill);


	// this even necessary? No idea
	G4int snum = step->GetNumberOfSecondariesInCurrentStep();
	if(snum == 0){return false;} //want this?
	
	const std::vector<const G4Track*>* secs = step->GetSecondaryInCurrentStep();

	LightDetectorHit** LDH = (LightDetectorHit**)malloc(snum*sizeof(LightDetectorHit*)); 

	//TODO: add primary vertices to event for secondaries
	//G4Event* sevent = G4RunManager::GetRunManager()->GetCurrentEvent();
	/////////////////
	// Was supposed to be for secondaries but unnecessary, I think
	/////////////////
	for(G4int i = 0; i < snum; i++){
	
		const G4Step* sstep = (*secs)[i]->GetStep();
		if(sstep == NULL){continue;}
	
		// track = (*secs)[i];
		populateHit(sstep, pdef,
			edep, stepLength, energy,
			pos, momdir,
			pdg, copynum, eventid, trackid, parentid,
			x,y,z,px,py,pz,
			partname,procname);


		//if(edep == 0.){continue;}  
		if (stepLength == 0.) {continue;}

		LDH[i] = new LightDetectorHit();

		edep = energy;
 
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

		partname = step->GetTrack()->GetParticleDefinition()->GetParticleName();
		
		LDH[i]->SetParticleName(partname);
		LDH[i]->SetProcessName(procname);

		// convert units m -> mm
		//LDH[i]->SetiPos(G4ThreeVector(stof(x)*1000,stof(y)*1000,stof(z)*1000));

		//LDH[i]->SetiMom(G4ThreeVector(stof(px),stof(py),stof(pz)));

		fHitsCollection->insert(LDH[i]);
	
		sstep->GetTrack()->SetTrackStatus(fStopAndKill);
	}


  return true;
}


void ELightDetectorSD::EndOfEvent(G4HCofThisEvent*)
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




