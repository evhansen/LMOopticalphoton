
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleTable.hh"

#include <random>
#include <chrono>


class G4Event;
class PrimaryGeneratorMessenger;
class DetectorConstruction;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  PrimaryGeneratorAction(DetectorConstruction*);
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleSource() { return fParticleSource; };

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);
  void SetRandomDirection(G4bool val = true);
  void SetRandomInLMO1(G4bool val = true);
  G4bool GetPolarized() { return fPolarized; };
  G4double GetPolarization() { return fPolarization; }

 private:

  char prng; // if 'c' use C PRNG else use G4 PRNG
  char misc;

  G4ParticleGun* fParticleSource;


  // *********
  // Define functions for computing the initial 
  // position and momenta of particles 
  // (qualities of the source)
  // *********
  
  //note: tot is the total number of iterations of the process 
  // to undergo and num is the number of particles generated 
  // in one iteration of the process; TODO: get rid of num 

  void RectanglePlateSource(G4int tot,G4int num,G4double PSHalfSizes[],G4double PSPos[],void (*Process)(float[],float[],G4Event*,G4ParticleGun*,char),G4Event* anEvent,char swch = 0);
  
  void CircularPlateSource(G4int tot,G4int num,G4double Radius,G4double PSPos[],void (*Process)(float[],float[],G4Event*,G4ParticleGun*,char),G4Event* anEvent,char swch = 0);
 
  // *********
  // Define the constituents of a single 
  // iteration of the processes
  // *********
 
  static void OPs(float ConstituentPos[],float ConstituentDir[],G4Event* anEvent,G4ParticleGun*PG,char swch);

  static void Cobalt60(float ConstituentPos[],float ConstituentDir[],G4Event* anEvent,G4ParticleGun*PG,char swch);
  
  static void CosmogenicMuons(float ConstituentPos[],float ConstituentDir[],G4Event* anEvent,G4ParticleGun*PG,char swch);
 

  // Misc
 
  PrimaryGeneratorMessenger* fGunMessenger;
  DetectorConstruction* fDetector;
  G4bool fRandomDirection;
  G4bool fRandomInLMO1;
  G4bool fPolarized;
  G4double fPolarization;

};





#endif /*PrimaryGeneratorAction_h*/
