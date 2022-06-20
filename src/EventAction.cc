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
// TODO: init event things here
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
	//G4int eventID = event->GetEventID();
		
	for(int i=0; i<nofHits; i++){
		LightDetectorHit* absoHit = (*absoHC)[i];
		
		auto analysisManager = G4AnalysisManager::Instance();
	
		
		analysisManager->FillNtupleDColumn(0, absoHit->GetEdep());
		
		
		analysisManager->FillNtupleIColumn(1, absoHit->GetPhysVolNum());
	
		analysisManager->FillNtupleFColumn(2, event->GetPrimaryVertex(0)->GetPosition().getX());
		analysisManager->FillNtupleFColumn(3, event->GetPrimaryVertex(0)->GetPosition().getY());
		analysisManager->FillNtupleFColumn(4, event->GetPrimaryVertex(0)->GetPosition().getZ());

		analysisManager->FillNtupleFColumn(5, absoHit->GetxfPos());
		analysisManager->FillNtupleFColumn(6, absoHit->GetyfPos());
		analysisManager->FillNtupleFColumn(7, absoHit->GetzfPos());
		
		
		analysisManager->FillNtupleIColumn(8, event->GetPrimaryVertex(0)->GetPrimary()->GetPDGcode());
	
		
		analysisManager->AddNtupleRow();
	}
}
