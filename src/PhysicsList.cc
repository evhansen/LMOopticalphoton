#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
//#include "G4EmStandardPhysicsSS.hh"
#include "G4MuonicAtomDecayPhysics.hh"
#include "G4PhotoElectricEffect.hh"

//#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTableIterator.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"






PhysicsList::PhysicsList() 
: G4VModularPhysicsList() {

  ver = 1;
  SetVerboseLevel(ver);

  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;
  G4cout << "<<< Geant4 Physics List simulation engine: FTFP_BERT"<<G4endl;
  G4cout <<G4endl;



ReplacePhysics(new G4EmStandardPhysics_option4);
  G4OpticalPhysics* oP = new G4OpticalPhysics();	 
  oP->Configure(kScintillation,true);
  oP->Configure(kCerenkov,false);

RegisterPhysics(oP);



/*
  //defaultCutValue = 0.7*CLHEP::mm;  
  defaultCutValue = 0.7*CLHEP::mm;  


  // EM Physics
  //RegisterPhysics( new G4EmStandardPhysics(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays 
  RegisterPhysics( new G4DecayPhysics(ver) );

  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver) );

  // Hadron Physics
  RegisterPhysics(  new G4HadronPhysicsFTFP_BERT(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));
  
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));

  //RegisterPhysics(new G4MuonicAtomDecayPhysics());  

  //ReplacePhysics(new G4EmStandardPhysics_option4());
  ReplacePhysics(new G4EmLivermorePhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  G4OpticalPhysics* oP = new G4OpticalPhysics();	 
  oP->Configure(kScintillation,true);
  oP->Configure(kCerenkov,false);

  auto oPa = G4OpticalParameters::Instance();
  oPa->SetWLSTimeProfile("delta");

  // new options broke scintillation?
  oPa->SetScintYieldFactor(1.0);
  oPa->SetScintExcitationRatio(0.0); // 1.0?
  oPa->SetScintEnhancedTimeConstants(true);
  oPa->SetScintTrackSecondariesFirst(true);

  //oPa->SetScintTrackInfo(false);
  oPa->SetScintStackPhotons(true);
  oPa->SetScintFiniteRiseTime(false);
  oPa->SetScintByParticleType(false);
 
  oPa->SetCerenkovMaxPhotonsPerStep(100);
  oPa->SetCerenkovMaxBetaChange(10.0);
  oPa->SetCerenkovTrackSecondariesFirst(true);
  
  //oPa->SetCerenkovStackPhotons(true);

  RegisterPhysics(oP);

  //G4StepLimiterPhysics* sL = new G4StepLimiterPhysics();
  //sL->SetApplyToAll(true);
  //RegisterPhysics(sL);
*/
}


PhysicsList::~PhysicsList(){ 
}

/*
void PhysicsList::ConstructParticle(){
}

void PhysicsList::ConstructProcess(){
}
*/

void PhysicsList::SetCuts(){
  G4VUserPhysicsList::SetCuts();
}  
