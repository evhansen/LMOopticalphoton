#include "RunAction.hh"
#include "Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "G4OpBoundaryProcess.hh"

RunAction::RunAction()
: G4UserRunAction()
{
	G4RunManager::GetRunManager()->SetPrintProgress(1);
	auto analysisManager = G4AnalysisManager::Instance();
	
	G4cout << "Using " << analysisManager->GetType() << G4endl;
	
	analysisManager->SetVerboseLevel(1);
	
	// G4Exception: merging ntuples not applicable in sequential application
	//analysisManager->SetNtupleMerging(true);

	//https://geant4.web.cern.ch/sites/default/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html

	analysisManager->CreateNtuple("LDHit", "Edep, SD, and Position");
	
	analysisManager->CreateNtupleFColumn("iEnergyMeV"); 
	
	analysisManager->CreateNtupleDColumn("EabsMeV");
 
	analysisManager->CreateNtupleIColumn("PhysVolNum"); 

	analysisManager->CreateNtupleFColumn("xiPosmm"); 
	analysisManager->CreateNtupleFColumn("yiPosmm"); 
	analysisManager->CreateNtupleFColumn("ziPosmm"); 

	// hit class
	analysisManager->CreateNtupleFColumn("xfPoshmm"); 
	analysisManager->CreateNtupleFColumn("yfPoshmm"); 
	analysisManager->CreateNtupleFColumn("zfPoshmm"); 
	
	// from event primary vertex
	analysisManager->CreateNtupleFColumn("xiPosPVmm"); 
	analysisManager->CreateNtupleFColumn("yiPosPVmm"); 
	analysisManager->CreateNtupleFColumn("ziPosPVmm"); 

	analysisManager->CreateNtupleIColumn("ParticleID"); 
	analysisManager->CreateNtupleIColumn("ProcessID"); 

	analysisManager->CreateNtupleFColumn("xiMom"); 
	analysisManager->CreateNtupleFColumn("yiMom"); 
	analysisManager->CreateNtupleFColumn("ziMom"); 

	analysisManager->CreateNtupleFColumn("xfMom"); 
	analysisManager->CreateNtupleFColumn("yfMom"); 
	analysisManager->CreateNtupleFColumn("zfMom"); 

	analysisManager->CreateNtupleIColumn("EventID"); 
	analysisManager->CreateNtupleIColumn("TrackID"); 
	analysisManager->CreateNtupleIColumn("ParentID"); 
	

	analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
	delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	auto analysisManager = G4AnalysisManager::Instance();
	G4String fileName = "";
	analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();
}
