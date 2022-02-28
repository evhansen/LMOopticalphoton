#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector)
: G4VUserActionInitialization(),fDetector(detector)
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
	SetUserAction(new RunAction);
}

void ActionInitialization::Build() const
{
	SetUserAction(new PrimaryGeneratorAction(fDetector));
	SetUserAction(new EventAction);
	SetUserAction(new RunAction);
}
