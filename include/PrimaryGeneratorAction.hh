#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class PrimaryGeneratorMessenger;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  PrimaryGeneratorAction(DetectorConstruction*);
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun() { return fParticleGun; };

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);
  void SetRandomDirection(G4bool val = true);
  void SetRandomInLMO1(G4bool val = true);
  G4bool GetPolarized() { return fPolarized; };
  G4double GetPolarization() { return fPolarization; }

 private:
  G4ParticleGun* fParticleGun;
  PrimaryGeneratorMessenger* fGunMessenger;
  DetectorConstruction* fDetector;
  G4bool fRandomDirection;
  G4bool fRandomInLMO1;
  G4bool fPolarized;
  G4double fPolarization;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/


// /// \file PrimaryGeneratorAction.hh
// /// \brief Definition of the PrimaryGeneratorAction class
//
// #ifndef PrimaryGeneratorAction_h
// #define PrimaryGeneratorAction_h 1
//
// #include "G4VUserPrimaryGeneratorAction.hh"
// #include "globals.hh"
//
// class G4ParticleGun;
// class G4Event;
//
// /// The primary generator action class with particle gum.
// ///
// /// It defines a single particle which hits the LightDetector
// /// perpendicular to the input face. The type of the particle
// /// can be changed via the G4 build-in commands of G4ParticleGun class
// /// (see the macros provided with this example).
//
// class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
// {
// public:
//   PrimaryGeneratorAction();
//   virtual ~PrimaryGeneratorAction();
//
//   virtual void GeneratePrimaries(G4Event* event);
//
//   // set methods
//   void SetRandomFlag(G4bool value);
//
// private:
//   G4ParticleGun*  fParticleGun; // G4 particle gun
// };
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// #endif
