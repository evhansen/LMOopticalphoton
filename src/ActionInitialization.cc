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
	
	RunAction* rA = new RunAction;
	SetUserAction(rA);

	// TODO: pass runaction to eventaction
	// in order to fix edep 	
	SetUserAction(new EventAction);
}
