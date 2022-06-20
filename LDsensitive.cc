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
/// \file example.cc
/// \brief Main program of the  example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"


#include "G4PhysListFactory.hh"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "FTFP_BERT.hh"

#include "G4String.hh"
#include "G4Types.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4StepLimiterPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " example [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  // Evaluate arguments
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Detect interactive mode (if no macro provided) and define UI session
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
#endif


  // TODO: fix physics registration after particles are dandy ..................

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
	 
  physicsList->RegisterPhysics(new G4OpticalPhysics());
  physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4EmExtraPhysics());
  physicsList->RegisterPhysics(new G4DecayPhysics());
  physicsList->RegisterPhysics(new G4HadronElasticPhysics());
  physicsList->RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  physicsList->RegisterPhysics(new G4StoppingPhysics());
  physicsList->RegisterPhysics(new G4IonPhysics());
  physicsList->RegisterPhysics(new G4NeutronTrackingCut());
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  
  runManager->SetUserInitialization(physicsList);

  
  // Set mandatory initialization classes
  
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  auto actionInitialization = new ActionInitialization(detConstruction);
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  
  G4VisManager* visManager = new G4VisExecutive(5);
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  
  if ( macro.size() ) { // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  
  } else  { // interactive mode : define UI session
    
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    
    UImanager->ApplyCommand("/control/execute LDsensitive.mac"); // yes, I am lazy

    ui->SessionStart();

    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....