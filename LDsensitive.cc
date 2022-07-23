/// \file example.cc
/// \brief Main program of the  example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"


#include "G4PhysListFactory.hh"

//#include "G4OpticalProcessIndex.hh"
#include "G4OpticalPhysics.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option4.hh"
//#include "G4EmStandardPhysicsSS.hh"
#include "G4MuonicAtomDecayPhysics.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "FTFP_BERT.hh"
//#include "QGSP_BERT.hh"

#include "G4String.hh"
#include "G4Types.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

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


  G4VModularPhysicsList* physicsList = new FTFP_BERT;


  G4OpticalPhysics* oP = new G4OpticalPhysics();	 
  oP->Configure(kScintillation,true);
  oP->Configure(kCerenkov,false);


  auto oPa = G4OpticalParameters::Instance();


  oPa->SetWLSTimeProfile("delta");

  // new options broke scintillation?
  oPa->SetScintYieldFactor(1.0);
  oPa->SetScintExcitationRatio(0.0); // 1.0?
  oPa->SetScintTrackSecondariesFirst(true);
  //oPa->SetScintTrackInfo(true);
  //oPa->SetScintStackPhotons(true);

  //oPa->SetScintFiniteRiseTime(false);
  oPa->SetScintEnhancedTimeConstants(false);
  //oPa->SetScintByParticleType(false);

 
  oPa->SetCerenkovMaxPhotonsPerStep(100);
  oPa->SetCerenkovMaxBetaChange(10.0);
  oPa->SetCerenkovTrackSecondariesFirst(true);
  //oPa->SetCerenkovStackPhotons(true);


  physicsList->RegisterPhysics(oP);


  //physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics());  
  //physicsList->RegisterPhysics(new G4MuonicAtomDecayPhysics());  
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());


  G4StepLimiterPhysics* sL = new G4StepLimiterPhysics();
  sL->SetApplyToAll(true);
  physicsList->RegisterPhysics(sL);


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
    UImanager->ApplyCommand("/control/execute "+macro);
  
  } else  { // interactive mode : define UI session
    
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    
    //UImanager->ApplyCommand("/control/execute LDsensitive.mac"); // yes, I am lazy

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
