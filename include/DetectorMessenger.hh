#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger : public G4UImessenger
{
 public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

 private:
  DetectorConstruction* fDetector;

  G4UIdirectory* fOpticalDir;

  // the LMO
  G4UIcmdWithAString* fLMOMatPropVectorCmd;
  G4UIcmdWithAString* fLMOMatPropConstCmd;
  G4UIcmdWithAString* fLMOMaterialCmd;

  // the LD
  G4UIcmdWithAString* fLDMatPropVectorCmd;
  G4UIcmdWithAString* fLDMatPropConstCmd;
  G4UIcmdWithAString* fLDMaterialCmd;

  // the world
  G4UIcmdWithAString* fWorldMatPropVectorCmd;
  G4UIcmdWithAString* fWorldMatPropConstCmd;
  G4UIcmdWithAString* fWorldMaterialCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
