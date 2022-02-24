#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class PrimaryGeneratorMessenger;
class DetectorConstruction;

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

  char pg; // if 1, using photon gun
  char pLD5; // if 1, pointing at LD5
  char zinc; // if 1, 0 incidence angle

  G4ParticleGun* fParticleGun;
  PrimaryGeneratorMessenger* fGunMessenger;
  DetectorConstruction* fDetector;
  G4bool fRandomDirection;
  G4bool fRandomInLMO1;
  G4bool fPolarized;
  G4double fPolarization;
};

#endif /*PrimaryGeneratorAction_h*/

