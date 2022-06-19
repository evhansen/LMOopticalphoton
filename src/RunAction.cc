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
	analysisManager->SetNtupleMerging(true);
	
	analysisManager->CreateNtuple("LDHit", "Edep, SD, and Position");
	
	
	analysisManager->CreateNtupleDColumn("Eabs"); 
	analysisManager->CreateNtupleIColumn("PhysVolNum"); 

	analysisManager->CreateNtupleFColumn("xiPos"); 
	analysisManager->CreateNtupleFColumn("yiPos"); 
	analysisManager->CreateNtupleFColumn("ziPos"); 

	analysisManager->CreateNtupleFColumn("xfPos"); 
	analysisManager->CreateNtupleFColumn("yfPos"); 
	analysisManager->CreateNtupleFColumn("zfPos"); 
	
	analysisManager->CreateNtupleIColumn("PDGID"); 

	//TODO: add parent particle for secondaries


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
