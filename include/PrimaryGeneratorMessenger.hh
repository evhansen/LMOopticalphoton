//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger : public G4UImessenger
{
 public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  virtual ~PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

 private:
  PrimaryGeneratorAction* fPrimaryAction;
  G4UIdirectory* fGunDir;
  G4UIcmdWithADoubleAndUnit* fPolarCmd;
  G4UIcmdWithABool* fRandomDirectionCmd;
  G4UIcmdWithABool* fRandomInLMO1Cmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
