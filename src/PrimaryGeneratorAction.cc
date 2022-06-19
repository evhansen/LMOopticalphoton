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


#include <random>
#include <chrono>


#include "DefSrcs.hh"

static void appendP(std::string filename, std::string line);
static char readExtFile(const char* filename);

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction() , fParticleSource(0), fDetector(det)
{

	prng = readExtFile("./options/prng");


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
	// Co-60 source, 37000/ s * Area
	//0.025 m diametre source => 18.162338 / s
	
	// Muon source, 4*10^(-7)/ cm^(2)s
	// 0.0225 ^3 m crystal => 0.03^2 m area => 3^2 cm area => 3.6 * 10^(-6) /s
	
	// => 10^6 seconds per segment so that:
	// 18162338 betas/ segment BUT upper half only => 9081169 betas/ segment BUT *2 for muons => 18162338
	// 3.6 muons / segment => *2 => ~7
		
	//99.88% of the radiation will be a 0.31 MeV e-, 1.17 MeV g, and 1.33 MeV g.

	//The muons will be 0.2-1000 GeV
	//ratio of m+/m- ~ 1.268 +/- (0.008 + 0.0002 p) where p is momentum.


	

	// TODO: fix this .. damn units	
	G4double LMOHS = fDetector->GetLMOXSize(); // LMOHalfSizes[0], in metres ..
	
	//sizes in meters
	G4double SquarePSHalfSizes[2] = {0.0225,0.0225};
	G4double Radius = 0.025;
	G4double SquarePSPos[3] = {0,0,0.04+0.0225};
	G4double CirclePSPos[3] = {0,0,-0.04-0.0225};



	
	// source 1: Muons
	//The muons will be 0.2-1000 GeV
	//note: hardcoded vertical
	//RectanglePlateSource(7,1,SquarePSHalfSizes,SquarePSPos,CosmogenicMuons,anEvent);

	// source 2: Co-60
	//99.88% of the radiation will be a 0.31 MeV e-, 1.17 MeV g, and 1.33 MeV g.
	//note: hardcoded rand unit vecs
	//CircularPlateSource(18162338,3,Radius,CirclePSPos,&Cobalt60,anEvent);



	// Life sux so reducing time segment significantly ...
	static unsigned int count = 0;

	CircularPlateSource(37,3,Radius,CirclePSPos,Cobalt60,anEvent);
	count++;

	if(count >= 490874){
		
		RectanglePlateSource(7,1,SquarePSHalfSizes,SquarePSPos,CosmogenicMuons,anEvent);

		count = 0;
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


// TODO: figure out if these are needed

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







void appendP(std::string filename, std::string line){
	std::ofstream p;
	p.open(filename.c_str(),std::ios_base::app);
	p << line;
	p.close();
}


char readExtFile(const char* filename){
	FILE *p = fopen(filename,"r");
	if(p == NULL){return '0';}
	return fgetc(p);
}
