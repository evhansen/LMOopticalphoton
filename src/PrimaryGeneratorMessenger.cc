#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* Gun)
  : G4UImessenger()
  , fPrimaryAction(Gun)
{
  fGunDir = new G4UIdirectory("/LDsensitive/gun/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fPolarCmd =
    new G4UIcmdWithADoubleAndUnit("/LDsensitive/gun/optPhotonPolar", this);
  fPolarCmd->SetGuidance("Set linear polarization");
  fPolarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  fPolarCmd->SetParameterName("angle", true);
  fPolarCmd->SetUnitCategory("Angle");
  fPolarCmd->SetDefaultValue(-360.0);
  fPolarCmd->SetDefaultUnit("deg");
  fPolarCmd->AvailableForStates(G4State_Idle);

  fRandomDirectionCmd =
    new G4UIcmdWithABool("/LDsensitive/gun/randomDirection", this);
  fRandomDirectionCmd->AvailableForStates(G4State_Idle, G4State_PreInit);

  fRandomInLMO1Cmd =
    new G4UIcmdWithABool("/LDsensitive/gun/randomInLMO1", this);
  fRandomInLMO1Cmd->AvailableForStates(G4State_Idle, G4State_PreInit);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPolarCmd;
  delete fGunDir;
  delete fRandomDirectionCmd;
  delete fRandomInLMO1Cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if(command == fPolarCmd)
  {
    G4double angle = fPolarCmd->GetNewDoubleValue(newValue);
    if(angle == -360.0 * deg)
    {
      fPrimaryAction->SetOptPhotonPolar();
    }
    else
    {
      fPrimaryAction->SetOptPhotonPolar(angle);
    }
  }
  else if(command == fRandomDirectionCmd)
  {
    fPrimaryAction->SetRandomDirection(true);
  }
  else if(command == fRandomInLMO1Cmd)
  {
	  fPrimaryAction->SetRandomInLMO1(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
