//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file LightDetectorSD.cc
/// \brief Implementation of the LightDetectorSD class

#include "LightDetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LightDetectorSD::~LightDetectorSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
  // for (G4int i=0; i<fNofCells+1; i++ ) {
  //   fHitsCollection->insert(new LightDetectorHit());
  // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LightDetectorSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{
  // energy deposit
  
  G4double edep = step->GetTotalEnergyDeposit();
  // was auto edep
 
  if(edep == 0.){
	/*G4cout<<"************************************************"<<G4endl;
  	G4cout<<"********NO ENERGY HIT DETECTED************"<<G4endl;
  	G4cout<<"************************************************"<<G4endl;

      	const G4VTouchable* touchable 
		= (step->GetPreStepPoint()->GetTouchable());
  
  	G4cout<<"Copy number: "<<touchable->GetCopyNumber()<<G4endl;
	*/
	
	  return false;  // No edep so don't count as hit
  }


  // step length
  // G4double stepLength = 0.;
  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
  //   stepLength = step->GetStepLength();
  // }
  // if (stepLength == 0. )
  //   return false; // No step length so don't count as hit

  G4StepPoint* thePrePoint = step->GetPreStepPoint();
  G4StepPoint* thePostPoint = step->GetPostStepPoint();

  const G4VTouchable* touchable = (step->GetPreStepPoint()->GetTouchable());
  // was auto touchable

  G4VPhysicalVolume* thePrePV = touchable->GetVolume();
  G4int copyNumber = touchable->GetCopyNumber();
  


    // G4cout<<"==> logical volume-> "<< thePrePV->GetLogicalVolume()->GetName() <<G4endl;
  // G4cout<<"==> which copy -> "<< touchable->GetCopyNumber() <<G4endl;

  G4ThreeVector pos = 
	  (thePrePoint->GetPosition() + thePostPoint->GetPosition()) / 2.;
  //pos /= 2.;


  // Get LightDetector cell id
  // auto layerNumber = touchable->GetReplicaNumber(1);
  //
  // // Get hit accounting data for this cell
  // auto hit = (*fHitsCollection)[layerNumber];
  // if ( ! hit ) {
  //   G4ExceptionDescription msg;
  //   msg << "Cannot access hit " << layerNumber;
  //   G4Exception("LightDetectorSD::ProcessHits()",
  //     "MyCode0004", FatalException, msg);
  // }

/* 
  G4cout<<"************************************************"<<G4endl;
  G4cout<<"************************************************"<<G4endl;
  G4cout<<"******************************\nHit Detected!! YAAAY\n******************************"<<G4endl;
  G4cout<<"************************************************"<<G4endl;
  G4cout<<"************************************************"<<G4endl;
  

  G4cout<<"Edep: "<< edep<<G4endl;
  G4cout<<"Copy number: "<< copyNumber<<G4endl;
  G4cout<<"Position: "<< pos<<G4endl;

  G4cout << "********* fHitsCollection->entries: "<<fHitsCollection->entries()
	  <<G4endl;
  for (G4int i=0;i<fHitsCollection->entries();i++){
  	G4cout<<"fHitsCollection entry: "<<i<<"\n Copy number: "
		<<(*fHitsCollection)[i]->GetPhysVolNum()<<G4endl;
  }
*/


  // Get hit for total accounting
  auto hitTotal
    = (*fHitsCollection)[fHitsCollection->entries()];

  LightDetectorHit* OPhit = new LightDetectorHit(thePrePV);

  OPhit->SetEdep(edep);
  OPhit->SetPhysVolNum(copyNumber);
  OPhit->SetPos(pos);
  
  //OPhit->SetiPos();

  fHitsCollection->insert(OPhit);

  // Add values
  // hit->Add(edep, stepLength);
  // hitTotal->Add(edep, stepLength);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
