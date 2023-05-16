#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include <sstream>
#include <iostream>



//TODO: add everything back, reimplement, etc.  without breaking things again

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
: G4UImessenger()
, fDetector(Det)
{
}

DetectorMessenger::~DetectorMessenger()
{
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
}
