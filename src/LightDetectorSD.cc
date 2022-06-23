#include "LightDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

LightDetectorSD::LightDetectorSD(
                            const G4String& name,
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

LightDetectorSD::~LightDetectorSD()
{
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
  
  G4double edep = step->GetTotalEnergyDeposit();
 
  if(edep == 0.){return false;}

  G4double stepLength = step->GetStepLength();
  
  if (stepLength == 0.) {return false;}

  G4StepPoint* thePrePoint = step->GetPreStepPoint();
  G4StepPoint* thePostPoint = step->GetPostStepPoint();

  G4ThreeVector pos = 
  	  (thePrePoint->GetPosition() + thePostPoint->GetPosition()) / 2.;

  G4Track* track = step->GetTrack(); // TODO: change engs?

  LightDetectorHit* OPhit = new LightDetectorHit();
 
  OPhit->SetEdep(edep);
  OPhit->SetPhysVolNum(thePrePoint->GetTouchable()->GetCopyNumber());
  OPhit->SetPos(pos);
  OPhit->SetTrackID(track->GetTrackID());
  OPhit->SetParentID(track->GetParentID());

  fHitsCollection->insert(OPhit);

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

