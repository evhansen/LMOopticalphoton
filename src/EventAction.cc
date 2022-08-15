#include "EventAction.hh"
#include "LightDetectorSD.hh"
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
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"


static int ParticleID(const G4ParticleDefinition*);


EventAction::EventAction()
: G4UserEventAction()
,fAbsHCID(-1)
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

	if ( fAbsHCID == -1 ) {
		fAbsHCID = G4SDManager::GetSDMpointer()->
			GetCollectionID("AbsorberHitsCollection");
	}
	
	LightDetectorHitsCollection* absoHC 
		= GetHitsCollection(fAbsHCID, event);
	
	G4int nofHits = absoHC->entries();

	for(int i=0; i<nofHits; i++){

		LightDetectorHit* absoHit = (*absoHC)[i];

		auto analysisManager = G4AnalysisManager::Instance();
	
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
	
		analysisManager->FillNtupleIColumn(12, ParticleID(absoHit->GetPDef()));

		analysisManager->FillNtupleFColumn(13, absoHit->GetxiMom());
		analysisManager->FillNtupleFColumn(14, absoHit->GetyiMom());
		analysisManager->FillNtupleFColumn(15, absoHit->GetziMom());

		analysisManager->FillNtupleFColumn(16, absoHit->GetxfMom());
		analysisManager->FillNtupleFColumn(17, absoHit->GetyfMom());
		analysisManager->FillNtupleFColumn(18, absoHit->GetzfMom());

		analysisManager->FillNtupleIColumn(19, absoHit->GetEventID());
		analysisManager->FillNtupleIColumn(20, absoHit->GetTrackID());
		analysisManager->FillNtupleIColumn(21, absoHit->GetParentID());
		
		analysisManager->AddNtupleRow();
	}
}





static int ParticleID(const G4ParticleDefinition* pd){


	if(pd == G4Electron::Definition()){
		return -11;
	} else if(pd == G4Positron::Definition()){
		return 11;
	} else if(pd == G4MuonMinus::Definition()){
		return -13;
	} else if(pd == G4MuonPlus::Definition()){
		return 13;
	} else if(pd ==G4Gamma::Definition()){
		return 22;
	} else if(pd ==G4OpticalPhoton::Definition()){
		return -22;
	} 

	return 0;		 
	
}

