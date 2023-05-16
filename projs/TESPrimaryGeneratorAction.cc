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


/////////////////
// Handle escaped OPs
/////////////////
void OPescape(G4double x, G4double y, G4double z, G4double px, G4double py, G4double pz, G4Event* anEvent, G4ParticleGun*fParticleSource){

	//G4cout << "OP escaped" << G4endl;

	fParticleSource->SetParticleEnergy(2.3*eV);
	fParticleSource->SetParticleDefinition(G4OpticalPhoton::Definition());
	fParticleSource->SetParticlePosition(G4ThreeVector(x*m,y*m,z*m));
	fParticleSource->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));


	if(fParticleSource->GetParticleDefinition() ==
		//G4OpticalPhoton::OpticalPhotonDefinition() ||
		G4OpticalPhoton::Definition())
	{
		if(fParticleSource->GetParticleDefinition() !=
		G4OpticalPhoton::OpticalPhotonDefinition())
		{
			G4ExceptionDescription ed;
			ed << "The particleGun is not an opticalphoton.";
			G4Exception("PrimaryGeneratorAction::SetOptPhotonPolar", "OpNovice2_004",JustWarning, ed);
			return;
		}

		G4bool fPolarized    = true;
		G4double fPolarization = 0.0 * deg;
		G4ThreeVector normal(1., 0., 0.);
		G4ThreeVector kphoton = fParticleSource->GetParticleMomentumDirection();
		G4ThreeVector product = normal.cross(kphoton);
		G4double modul2       = product * product;
		G4ThreeVector e_perpend(0., 0., 1.);
	
		if(modul2 > 0.){
		e_perpend = (1. / std::sqrt(modul2)) * product;
		}
		G4ThreeVector e_paralle = e_perpend.cross(kphoton);
		G4ThreeVector polar =
		std::cos(fPolarization) * e_paralle + std::sin(fPolarization) * e_perpend;
		fParticleSource->SetParticlePolarization(polar);
	}


	fParticleSource->GeneratePrimaryVertex(anEvent);
}


PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction() , fParticleSource(0), fDetector(det)
{
	// set PRNG //
	prng = readExtFile("./options/prng");
	// misc used for troubleshotting LMO runs //
	misc = readExtFile("./options/misc");

	fGunMessenger = new PrimaryGeneratorMessenger(this);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	fRandomDirection = false;
	fPolarized       = false;
	fPolarization    = 0.;
	
	fParticleSource = new G4ParticleGun(1);

	/////////////////
	// Default particle source settings
	/////////////////
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
	/////////////////
	// Set options
	/////////////////
	
	// set whether OPs escape from wires
	static char esc = readExtFile("./options/TES_escape");
	
	// set the type of optical fibre being used (horizontal or vertical)
	static char ofc = readExtFile("./options/TES_opticalfibre");
	
	// set the top gap value (default 10mm -- value is (10 + file value) in mm)
	float tgap = stof(readExtFileLine("./options/TES_topgap"));



	/////////////////
	// Setting the geometry of the fibres
	/////////////////
	
	G4double Mm = 0.001;

	tgap = tgap * Mm;
	
	G4double FRAMESIDE = (55.5/2)*Mm;
	G4double SIDELONG = (135/2)*Mm; // 127mm + 8mm top slab thick
	G4double TOPLONG = (102/2)*Mm;

	G4double SideHalfSizes[3] = {SIDELONG + tgap, FRAMESIDE, 2*Mm};
	G4double SlatHalfSizes[3] = {TOPLONG, FRAMESIDE, 3.5*Mm};
	G4double TopHalfSizes[3] = {TOPLONG, FRAMESIDE, 4*Mm};

	G4double SlatInit = 17*Mm - SlatHalfSizes[2];
	G4double SlatHeightDiff = 25*Mm;
	G4double FibreInit = (SlatHeightDiff/2)+SlatInit;
	G4double LastSlatPos = SlatInit+(4*SlatHeightDiff);
	G4double TopPos = SideHalfSizes[0]*2 - TopHalfSizes[2];

	// HORIZONTAL FIBRE VALUES //
	G4double HHS[3] = {0.102/2,0.000000005,0.000000005};
	G4double HPos1[3] = {0,0,FibreInit};
	G4double HPos2[3] = {0,0,FibreInit+SlatHeightDiff};
	G4double HPos3[3] = {0,0,FibreInit+(2*SlatHeightDiff)};
	G4double HPos4[3] = {0,0,FibreInit+(3*SlatHeightDiff)};
	G4double HPos5[3] = {0,0,((TopPos-LastSlatPos)/2)+LastSlatPos};
	
	// VERTICAL FIBRE VALUES //
	G4double VHS[3] = {0.000000005,0.000000005,0.135/2};
	G4double VPos1[3] = {0.0245,0.03,0.135/2};
	G4double VPos2[3] = {-0.0245,0.03,0.135/2};



	static int cnt = 1;

	// Horizontal
	if(ofc == 'h'){
		//horiz
		switch(cnt){
			case 1:
			RectanglePlateSource(1,1,HHS,HPos1,OPs,anEvent); // define source; see documentation
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(102/2,HPos1[1],HPos1[2],-0.9,-0.,0,anEvent,fParticleSource);return;} // escape
			break;

			case 2:
			RectanglePlateSource(1,1,HHS,HPos2,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(102/2,HPos2[1],HPos2[2],-0.9,-0.,0,anEvent,fParticleSource);return;}
			break;

			case 3:
			RectanglePlateSource(1,1,HHS,HPos3,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(102/2,HPos3[1],HPos3[2],-0.9,-0.,0,anEvent,fParticleSource);return;}
			break;

			case 4:
			RectanglePlateSource(1,1,HHS,HPos4,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(102/2,HPos4[1],HPos4[2],-0.9,-0.,0,anEvent,fParticleSource);return;}
			break;

			case 5:
			RectanglePlateSource(1,1,HHS,HPos5,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(102/2,HPos5[1],HPos5[2],-0.9,-0.,0,anEvent,fParticleSource);return;}
			break;

			default:
			cnt = 0;
			break;
		}

		cnt++;
		
		if(cnt < 1 || cnt > 5){
			cnt = 1;
		}



	//Vertical
	} else if (ofc == 'v') {
		//in front
		switch(cnt){
			case 1:
			RectanglePlateSource(1,1,VHS,VPos1,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(VPos1[0],VPos1[1],0,0,-0.05,-0.9,anEvent,fParticleSource);return;}
			break;

			case 2:
			RectanglePlateSource(1,1,VHS,VPos2,OPs,anEvent);
			if(esc == 'a' && ((double)rand()/(double)RAND_MAX) < 0.3){OPescape(VPos2[0],VPos2[1],0,0,-0.05,-0.9,anEvent,fParticleSource);return;}
			break;

			default:
			cnt = 0;
			break;
		}

		cnt++;
		
		if(cnt < 1 || cnt > 2){
			cnt = 1;
		}



// Add different optical fibre geometry
//	} else { 
//		//along walls

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



