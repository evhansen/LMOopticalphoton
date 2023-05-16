#include "EventAction.hh"

//#include "LightDetectorSD.hh"
//#include "ELightDetectorSD.hh"
#include "LightDetectorHit.hh"
#include "Analysis.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Alpha.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"


static int ParticleID(const G4ParticleDefinition* pd, G4String pn);
static int ProcessID(const G4String name);



EventAction::EventAction()
: G4UserEventAction()
,fAbsHCID(-1)
,fHCID(-1)
{}

EventAction::~EventAction()
{}

LightDetectorHitsCollection*
EventAction::GetHitsCollection(G4int hcID,
const G4Event* event) const
{

	auto hitsCollection
	= static_cast<LightDetectorHitsCollection*>(
	event->GetHCofThisEvent()->GetHC(hcID));
		
	if ( ! hitsCollection ) {
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID " << hcID;
		G4Exception("EventAction::GetHitsCollection()",
		"MyCode0003", FatalException, msg);
	}
	
	return hitsCollection;
}

void EventAction::PrintEventStatistics(
G4double absoEdep//, G4double absoTrackLength
) const
{
	G4cout
	<< "   Absorber: total energy: "
	<< std::setw(7) << G4BestUnit(absoEdep, "Energy")
	<< G4endl;
}

void EventAction::BeginOfEventAction(const G4Event*)
{
}

void EventAction::EndOfEventAction(const G4Event* event)
{

	if ( fHCID == -1 ) {
		fHCID = G4SDManager::GetSDMpointer()->
			GetCollectionID("HitsCollection");
	}

	if ( fAbsHCID == -1 ) {
		fAbsHCID = G4SDManager::GetSDMpointer()->
			GetCollectionID("AbsorberHitsCollection");
	}

	/////////////////
	// From the regular sensitive detectors
	/////////////////
	LightDetectorHitsCollection* soHC 
		= GetHitsCollection(fHCID, event);
	
	/////////////////
	// From the absorber sensitive detectors
	/////////////////
	LightDetectorHitsCollection* absoHC 
		= GetHitsCollection(fAbsHCID, event);

	auto analysisManager = G4AnalysisManager::Instance();

	
	

	//for(int i = 0; i < soHC->entries(); i++){
	//	LightDetectorHit* absoHit = (*soHC)[i];

	for(int i = 0; i < absoHC->entries() + soHC->entries(); i++){

		LightDetectorHit* absoHit;
	
		if(i<absoHC->entries()){
			absoHit = (*absoHC)[i];
		}else{
			absoHit = (*soHC)[i-absoHC->entries()];
		}

		// need to ditch things that have ridiculous energies, 
		// so particles from weird procs that don't make sense
		if(ParticleID(absoHit->GetPDef(),absoHit->GetParticleName()) == 0){

			std::string write = "Testing: " + std::to_string(absoHit->GetEnergy()) +
				std::to_string(absoHit->GetEdep()) +
				std::to_string(absoHit->GetPhysVolNum()) +
				std::to_string(absoHit->GetxiPos()) +
				std::to_string(absoHit->GetyiPos()) +
				std::to_string(absoHit->GetziPos()) +
				std::to_string(absoHit->GetxfPos()) +
				std::to_string(absoHit->GetyfPos()) +
				std::to_string(absoHit->GetzfPos()) +
				std::to_string(event->GetPrimaryVertex(0)->GetPosition().getX()) +
				std::to_string(event->GetPrimaryVertex(0)->GetPosition().getY()) +
				std::to_string(event->GetPrimaryVertex(0)->GetPosition().getZ()) +
				std::to_string(ParticleID(absoHit->GetPDef(),absoHit->GetParticleName())) +
				std::to_string(ProcessID(absoHit->GetProcessName())) +
				std::to_string(absoHit->GetxiMom()) +
				std::to_string(absoHit->GetyiMom()) +
				std::to_string(absoHit->GetziMom()) +
				std::to_string(absoHit->GetxfMom()) +
				std::to_string(absoHit->GetyfMom()) +
				std::to_string(absoHit->GetzfMom()) +
				std::to_string(absoHit->GetEventID()) +
				std::to_string(absoHit->GetTrackID()) +
				std::to_string(absoHit->GetParentID()) + "\n";

				G4cout << write << G4endl;
		}	

		analysisManager->FillNtupleFColumn(0, absoHit->GetEnergy());
	
		analysisManager->FillNtupleDColumn(1, absoHit->GetEdep());
	
		analysisManager->FillNtupleIColumn(2, absoHit->GetPhysVolNum());

		analysisManager->FillNtupleFColumn(3, absoHit->GetxiPos());
		analysisManager->FillNtupleFColumn(4, absoHit->GetyiPos());
		analysisManager->FillNtupleFColumn(5, absoHit->GetziPos());

		analysisManager->FillNtupleFColumn(6, absoHit->GetxfPos());
		analysisManager->FillNtupleFColumn(7, absoHit->GetyfPos());
		analysisManager->FillNtupleFColumn(8, absoHit->GetzfPos());	

		analysisManager->FillNtupleFColumn(9, event->GetPrimaryVertex(0)->GetPosition().getX());
		analysisManager->FillNtupleFColumn(10, event->GetPrimaryVertex(0)->GetPosition().getY());
		analysisManager->FillNtupleFColumn(11, event->GetPrimaryVertex(0)->GetPosition().getZ());
	
		analysisManager->FillNtupleIColumn(12, ParticleID(absoHit->GetPDef(),absoHit->GetParticleName()));
		analysisManager->FillNtupleIColumn(13, ProcessID(absoHit->GetProcessName()));
		
		analysisManager->FillNtupleFColumn(14, absoHit->GetxiMom());
		analysisManager->FillNtupleFColumn(15, absoHit->GetyiMom());
		analysisManager->FillNtupleFColumn(16, absoHit->GetziMom());

		analysisManager->FillNtupleFColumn(17, absoHit->GetxfMom());
		analysisManager->FillNtupleFColumn(18, absoHit->GetyfMom());
		analysisManager->FillNtupleFColumn(19, absoHit->GetzfMom());

		analysisManager->FillNtupleIColumn(20, absoHit->GetEventID());
		analysisManager->FillNtupleIColumn(21, absoHit->GetTrackID());
		analysisManager->FillNtupleIColumn(22, absoHit->GetParentID());


		//G4cout << "TESTING: " << absoHit->GetParticleName() << "\n" << absoHit->GetProcessName() << "\n" << G4endl; 
		
		analysisManager->AddNtupleRow();
	}
}





static int ParticleID(const G4ParticleDefinition* pd, G4String pn){
	if(pd == G4Electron::Definition()){
		return -11;
	} else if(pd == G4Alpha::Definition()){
		return 12;
	} else if(pd == G4Positron::Definition()){
		return 11;
	} else if(pd == G4MuonMinus::Definition()){
		return -13;
	} else if(pd == G4MuonPlus::Definition()){
		return 13;
	} else if(pd == G4Gamma::Definition()){
		return 22;
	} else if(pd == G4OpticalPhoton::Definition()){
		return -22;
	}

	if(pn.empty()){return -1;}
	
	if(pn == "e-"){
		return -11;
	} else if(pn == "e+"){
		return 11;
	} else if(pn == "gamma"){
		return 22;
	} else if(pn == "opticalphoton"){
		return -22;
	}

	//G4cout << "WARNING: " << pd->GetParticleName() << "\n" << G4endl; 
	G4cout << "WARNING: " << pn << "\n" << G4endl; 
	
	G4cout << "WARNING: ParticleID.\n" << G4endl; 
	
	return 0;		 
}

static int ProcessID(const G4String name){
	if(name == "msc"){ // e-
		return 1;
	} else if(name == "eIoni"){ // e-
		return 2;
	} else if(name == "eBrem"){ // e-
		return 3;
	} else if(name == "OpAbsorption"){ //OP
		return 11;
	} else if(name == "compt"){ // gamma
		return 21;
	} else if(name == "phot"){ // gamma
		return 22;
	} else if(name == "Transportation") {
		return 31;
	}

	if(name.empty()){return -1;}

	G4cout << "WARNING: " << name << "\n" << G4endl; 
	G4cout << "WARNING: ProcessID.\n" << G4endl; 

	return 0;		 
}

