#include "RunAction.hh"
#include "Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpBoundaryProcess.hh"

RunAction::RunAction()
: G4UserRunAction()
{
	G4RunManager::GetRunManager()->SetPrintProgress(1);
	auto analysisManager = G4AnalysisManager::Instance();
	
	G4cout << "Using " << analysisManager->GetType() << G4endl;
	
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true);
	
	analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 5*eV);
	analysisManager->CreateH1("PhysVolNum","Which SD is hit", 30, 0., 30);
	
	analysisManager->CreateNtuple("LDHit", "Edep, SD, and Position");
	
	
	analysisManager->CreateNtupleDColumn("Eabs"); 
	analysisManager->CreateNtupleIColumn("PhysVolNum"); 

	analysisManager->CreateNtupleFColumn("xiPos"); 
	analysisManager->CreateNtupleFColumn("yiPos"); 
	analysisManager->CreateNtupleFColumn("ziPos"); 

	analysisManager->CreateNtupleFColumn("xfPos"); 
	analysisManager->CreateNtupleFColumn("yfPos"); 
	analysisManager->CreateNtupleFColumn("zfPos"); 


	// I suspect that D, I, F, etc. are double, integer, float, etc.

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
		
	if ( analysisManager->GetH1(0) ) {
		G4cout << G4endl << " ----> print histograms statistic ";
		
		if(isMaster) {
			G4cout << "for the entire run " << G4endl << G4endl;
		} else {
			G4cout << "for the local thread " << G4endl << G4endl;
		}

		G4cout << " EAbs : mean = "
		<< G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
		<< " count = "
		<< analysisManager->GetH1(0)->entries() << G4endl;
	}
	analysisManager->Write();
	analysisManager->CloseFile();
}
