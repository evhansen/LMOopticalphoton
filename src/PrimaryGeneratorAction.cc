#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"


#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"

#include "DefSrcs.hh"



PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction() , fParticleSource(0), fDetector(det)
{
	prng = readExtFile("./options/prng");
	misc = readExtFile("./options/misc");

	fGunMessenger = new PrimaryGeneratorMessenger(this);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	fRandomDirection = false;
	fPolarized       = false;
	fPolarization    = 0.;
	
	fParticleSource = new G4ParticleGun(1);

	fParticleSource->SetParticleEnergy(2.3*eV);
	fParticleSource->SetParticleDefinition(G4Electron::Definition());
	
	fParticleSource->SetParticleTime(0.0 * ns);
	fParticleSource->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
	fParticleSource->SetParticleMomentumDirection(G4ThreeVector(0., 0., 0.));



}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleSource;
	delete fGunMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	// TODO: fix this .. damn units	
	G4double LMOHS = fDetector->GetLMOXSize(); // LMOHalfSizes[0], in metres ..
	
	//sizes in meters
	G4double SquarePSHalfSizes[2] = {0.0225,0.0225};
	G4double Radius = 0.025;
	G4double SquarePSPos[3] = {0,0,0.04+0.0225};
	G4double CirclePSPos[3] = {0,0,-0.04-0.0225};

	if(misc == 'b' || misc == 'g' || misc == '1' || misc == '2'){
		CircularPlateSource(1,3,Radius,CirclePSPos,Cobalt60,anEvent,misc);
	} else {
		RectanglePlateSource(1,1,SquarePSHalfSizes,SquarePSPos,CosmogenicMuons,anEvent);
	}

	if(fParticleSource->GetParticleDefinition() ==
		//G4OpticalPhoton::OpticalPhotonDefinition() ||
		G4OpticalPhoton::Definition())
	{
		if(fPolarized)
		SetOptPhotonPolar(fPolarization);
		else
		SetOptPhotonPolar();
	}
}



void PrimaryGeneratorAction::SetOptPhotonPolar()
{
	//G4double angle = G4UniformRand() * 360.0 * deg;
	G4double angle = 0.0 * deg;
	SetOptPhotonPolar(angle);
}

void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
	if(fParticleSource->GetParticleDefinition() !=
	G4OpticalPhoton::OpticalPhotonDefinition())
	{
		G4ExceptionDescription ed;
		ed << "The particleGun is not an opticalphoton.";
		G4Exception("PrimaryGeneratorAction::SetOptPhotonPolar", "OpNovice2_004",
		JustWarning, ed);
		return;
	}

	fPolarized    = true;
	fPolarization = angle;
	G4ThreeVector normal(1., 0., 0.);
	G4ThreeVector kphoton = fParticleSource->GetParticleMomentumDirection();
	G4ThreeVector product = normal.cross(kphoton);
	G4double modul2       = product * product;
	G4ThreeVector e_perpend(0., 0., 1.);
	
	if(modul2 > 0.)
	e_perpend = (1. / std::sqrt(modul2)) * product;
	G4ThreeVector e_paralle = e_perpend.cross(kphoton);
	G4ThreeVector polar =
	std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
	fParticleSource->SetParticlePolarization(polar);
}

void PrimaryGeneratorAction::SetRandomDirection(G4bool val)
{
	fRandomDirection = val;
}

void PrimaryGeneratorAction::SetRandomInLMO1(G4bool val)
{
	fRandomInLMO1 = val;
}



