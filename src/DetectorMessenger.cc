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

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
: G4UImessenger()
, fDetector(Det)
{
	fOpticalDir = new G4UIdirectory("/LDsensitive/");
	fOpticalDir->SetGuidance("Parameters for optical simulation.");

	fLMOMatPropVectorCmd =
	new G4UIcmdWithAString("/LDsensitive/LMOProperty", this);

	fLMOMatPropVectorCmd->SetGuidance("Set material property vector for ");
	fLMOMatPropVectorCmd->SetGuidance("the LMO crystals.");
	fLMOMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLMOMatPropVectorCmd->SetToBeBroadcasted(false);

	fLMOMatPropConstCmd =
	new G4UIcmdWithAString("/LDsensitive/LMOConstProperty", this);

	fLMOMatPropConstCmd->SetGuidance("Set material constant property ");
	fLMOMatPropConstCmd->SetGuidance("for the LMO crystals.");
	fLMOMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLMOMatPropConstCmd->SetToBeBroadcasted(false);

	fLMOMaterialCmd = new G4UIcmdWithAString("/LDsensitive/LMOMaterial", this);

	fLMOMaterialCmd->SetGuidance("Set material of the LMO crystals.");
	fLMOMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLMOMaterialCmd->SetToBeBroadcasted(false);

	fLDMatPropVectorCmd =
	new G4UIcmdWithAString("/LDsensitive/LDProperty", this);

	fLDMatPropVectorCmd->SetGuidance("Set material property vector for ");
	fLDMatPropVectorCmd->SetGuidance("the light detectors.");
	fLDMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLDMatPropVectorCmd->SetToBeBroadcasted(false);

	fLMOMatPropConstCmd =
	new G4UIcmdWithAString("/LDsensitive/LDConstProperty", this);

	fLMOMatPropConstCmd->SetGuidance("Set material constant property ");
	fLMOMatPropConstCmd->SetGuidance("for the light detectors.");
	fLMOMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLMOMatPropConstCmd->SetToBeBroadcasted(false);

	fLMOMaterialCmd = new G4UIcmdWithAString("/LDsensitive/LDMaterial", this);

	fLMOMaterialCmd->SetGuidance("Set material of the light detectors.");
	fLMOMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLMOMaterialCmd->SetToBeBroadcasted(false);

	fWorldMatPropVectorCmd =
	new G4UIcmdWithAString("/LDsensitive/worldProperty", this);

	fWorldMatPropVectorCmd->SetGuidance("Set material property vector ");
	fWorldMatPropVectorCmd->SetGuidance("for the world.");
	fWorldMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fWorldMatPropVectorCmd->SetToBeBroadcasted(false);

	fWorldMatPropConstCmd =
	new G4UIcmdWithAString("/LDsensitive/worldConstProperty", this);

	fWorldMatPropConstCmd->SetGuidance("Set material constant property");
	fWorldMatPropConstCmd->SetGuidance(" for the world.");
	fWorldMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fWorldMatPropConstCmd->SetToBeBroadcasted(false);

	fWorldMaterialCmd = new G4UIcmdWithAString("/LDsensitive/worldMaterial", this);

	fWorldMaterialCmd->SetGuidance("Set material of world.");
	fWorldMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fWorldMaterialCmd->SetToBeBroadcasted(false);
}

DetectorMessenger::~DetectorMessenger()
{
	delete fOpticalDir;
	delete fLMOMatPropVectorCmd;
	delete fLMOMatPropConstCmd;
	delete fLMOMaterialCmd;
	delete fWorldMatPropVectorCmd;
	delete fWorldMatPropConstCmd;
	delete fWorldMaterialCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fLMOMatPropVectorCmd)
	{
		G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
		std::istringstream instring(newValue);
		G4String prop;
		instring >> prop;
		while(instring)
		{
			G4String tmp;
			instring >> tmp;
			if(tmp == "")
			{
				break;
			}
			G4double en = G4UIcommand::ConvertToDouble(tmp);
			instring >> tmp;
			G4double val = G4UIcommand::ConvertToDouble(tmp);
			mpv->InsertValues(en, val);
		}
		//fDetector->AddLMOMPV(prop, mpv);
	}
	else if(command == fLDMatPropVectorCmd)
	{
		G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
		std::istringstream instring(newValue);
		G4String prop;
		instring >> prop;
		while(instring)
		{
			G4String tmp;
			instring >> tmp;
			if(tmp == "")
			{
				break;
			}
			G4double en = G4UIcommand::ConvertToDouble(tmp);
			instring >> tmp;
			G4double val = G4UIcommand::ConvertToDouble(tmp);
			mpv->InsertValues(en, val);
		}
		//fDetector->AddLDMPV(prop, mpv);
	}
	else if(command == fWorldMatPropVectorCmd)
	{
		G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
		std::istringstream instring(newValue);
		G4String prop;
		instring >> prop;
		while(instring)
		{
			G4String tmp;
			instring >> tmp;
			if(tmp == "")
			{
				break;
			}
			G4double en = G4UIcommand::ConvertToDouble(tmp);
			instring >> tmp;
			G4double val = G4UIcommand::ConvertToDouble(tmp);
			mpv->InsertValues(en, val);
		}
		//fDetector->AddWorldMPV(prop, mpv);
	}
	else if(command == fLMOMatPropConstCmd)
	{
		std::istringstream instring(newValue);
		G4String prop;
		G4String tmp;
		instring >> prop;
		instring >> tmp;
		G4double val = G4UIcommand::ConvertToDouble(tmp);
		//fDetector->AddLMOMPC(prop, val);
	}
	else if(command == fLDMatPropConstCmd)
	{
		std::istringstream instring(newValue);
		G4String prop;
		G4String tmp;
		instring >> prop;
		instring >> tmp;
		G4double val = G4UIcommand::ConvertToDouble(tmp);
		//fDetector->AddLMOMPC(prop, val);
	}
	else if(command == fWorldMatPropConstCmd)
	{
		std::istringstream instring(newValue);
		G4String prop;
		G4String tmp;
		instring >> prop;
		instring >> tmp;
		G4double val = G4UIcommand::ConvertToDouble(tmp);
		//fDetector->AddWorldMPC(prop, val);
	}
	else if(command == fWorldMaterialCmd)
	{
		//fDetector->SetWorldMaterial(newValue);
	}
	else if(command == fLMOMaterialCmd)
	{
		//fDetector->SetLMOMaterial(newValue);
	}
	else if(command == fLDMaterialCmd)
	{
		//fDetector->SetLDMaterial(newValue);
	}
}
