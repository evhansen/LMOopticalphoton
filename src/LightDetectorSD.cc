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

  // Create hits
  // fNofCells for cells + one more for total sums
  for (G4int i=0; i<fNofCells+1; i++ ) {
  	fHitsCollection->insert(new LightDetectorHit());
  }
}


G4bool LightDetectorSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{
  // energy deposit
  
  G4double edep = step->GetTotalEnergyDeposit();
  // was auto edep
 
  if(edep == 0.){
	/*const G4VTouchable* touchable 
		= (step->GetPreStepPoint()->GetTouchable());
  
  	G4cout<<"Copy number: "<<touchable->GetCopyNumber()<<G4endl;
	*/
	
	  return false;  // No edep so don't count as hit
  }


  // step length
  /*G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
     stepLength = step->GetStepLength();
  }*/
  G4double stepLength = step->GetStepLength();
  
  if (stepLength == 0. ) {
	return false; // No step length so don't count as hit
  }

  G4StepPoint* thePrePoint = step->GetPreStepPoint();
  G4StepPoint* thePostPoint = step->GetPostStepPoint();

  const G4VTouchable* touchable = (step->GetPreStepPoint()->GetTouchable());

  G4VPhysicalVolume* thePrePV = touchable->GetVolume();
  G4int copyNumber = touchable->GetCopyNumber();
  

  G4ThreeVector pos = 
	  (thePrePoint->GetPosition() + thePostPoint->GetPosition()) / 2.;


  // Get LightDetector cell id
  // auto layerNumber = touchable->GetReplicaNumber(1);
  
  // Get hit accounting data for this cell
  // auto hit = (*fHitsCollection)[layerNumber];
  // if ( ! hit ) {
  //   G4ExceptionDescription msg;
  //   msg << "Cannot access hit " << layerNumber;
  //   G4Exception("LightDetectorSD::ProcessHits()",
  //     "MyCode0004", FatalException, msg);
  // }


  // Get hit for total accounting
  auto hitTotal
    = (*fHitsCollection)[fHitsCollection->entries()];

  LightDetectorHit* OPhit = new LightDetectorHit(thePrePV);

  OPhit->SetEdep(edep);
  OPhit->SetPhysVolNum(copyNumber);
  OPhit->SetPos(pos);

  fHitsCollection->insert(OPhit);

  // Add values
  //hit->Add(edep, stepLength);
  //hitTotal->Add(edep, stepLength);

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

