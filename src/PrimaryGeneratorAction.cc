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


static char readExtFile(const char* filename);

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction() , fParticleGun(0), fDetector(det)
{

	pg = readExtFile("./options/pg");
	pLD5 = readExtFile("./options/pld5");
	zinc = readExtFile("./options/zinc");
	
	G4int n_particle = 1;

	fParticleGun     = new G4ParticleGun(n_particle);

	fGunMessenger = new PrimaryGeneratorMessenger(this);

	fRandomDirection = false;
	fPolarized       = false;
	fPolarization    = 0.;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleTime(0.0 * ns);
	fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 0.));

	//fParticleGun->SetParticleEnergy(500.0 * keV);
	fParticleGun->SetParticleEnergy(2.6 * eV);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fGunMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4double vx, vy, vz, ux, uy, uz;

		if(pg != '1' && pg != '2') {
			G4double cosTheta = 2*G4UniformRand() - 1.;
			G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
			G4double phi = CLHEP::twopi*G4UniformRand();
			vx = sinTheta*std::cos(phi);
			vy = sinTheta*std::sin(phi);
			vz = cosTheta;
			fParticleGun->SetParticleMomentumDirection
				(G4ThreeVector(vx,vy,vz));

			// LMO half length
			G4double Rmax = fDetector->GetLMOXSize();
			
			// LMO length
			G4double R2 = Rmax*2.0;
		
			// G4UniformRand is in (0,1)
			// let LMO halflength = H
			// need affine ~U(0,1) => ~U(-H,H)
			// for x\in (0,1) => y = -H + x (H - (-H)) = -H + 2H x 
			//ux = G4UniformRand()*R2 - Rmax; 
			//uy = G4UniformRand()*R2 - Rmax; 
			//uz = G4UniformRand()*R2 - Rmax;

			float rnd1 = ((float)rand()/RAND_MAX);
			float rnd2 = ((float)rand()/RAND_MAX);
			float rnd3 = ((float)rand()/RAND_MAX);
			
			ux = rnd1*R2 - Rmax; 
			uy = rnd2*R2 - Rmax; 
			uz = rnd3*R2 - Rmax;

			fParticleGun->SetParticlePosition(G4ThreeVector(ux,uy,uz));
	
		} else if (pg == '2') {
				
			
			if(pLD5 == '1'){

				if(zinc == '1'){
					vx = -1.0;
					vy = 0.0;
					vz = 0.0;
				
					ux = 0.0*mm;
					uy = 0.0*mm;
					uz = 0.0*mm;
				} else {
				
					vx = -0.5;
					vy = -0.5 / std::tan(CLHEP::twopi/(9*2));
					vz = 0.0;
				
					ux = 0.0*mm;
					uy = 10.0*mm;
					uz = 0.0*mm;
				}

			
				fParticleGun->SetParticleMomentumDirection
					(G4ThreeVector(vx,vy,vz));
								
				fParticleGun->SetParticlePosition
					(G4ThreeVector(ux,uy,uz));


			} else {
				
	
				if(zinc == '1'){
					vx = 0.0;
					vy = 0.0;
					vz = 1.0;
				
					ux = 0.0*mm;
					uy = 0.0*mm;
					uz = 0.0*mm;
				} else {
				
					vy = 0.0;
					vx = -0.5 / std::tan(CLHEP::twopi/(9*2));
					vz = 0.5;
				
					uy = 0.0*mm;
					ux = 10.0*mm;
					uz = 0.0*mm;
				}

			
				fParticleGun->SetParticleMomentumDirection
					(G4ThreeVector(vx,vy,vz));
								
				fParticleGun->SetParticlePosition
					(G4ThreeVector(ux,uy,uz));


			}

		} else { // big boi for fibre
	
			vx = -1.0;
			vy = 0.0;
			//G4double vy = -1.0/std::tan(CLHEP::twopi/(9*2));
			vz = 0.0;
			fParticleGun->SetParticleMomentumDirection
				(G4ThreeVector(vx,vy,vz));

			ux = 33.0*mm; 
			uy = 0.0*mm;
			uz = 0.0*mm;
			fParticleGun->SetParticlePosition(G4ThreeVector(ux,uy,uz));
		}
	


		if(fParticleGun->GetParticleDefinition() ==
			G4OpticalPhoton::OpticalPhotonDefinition())
		{
			if(fPolarized)
			SetOptPhotonPolar(fPolarization);
			else
			SetOptPhotonPolar();
		}

	//G4cout //<< "Position and momentum generated\n"
		//<< "***********\nMomentum dir vxyz: " <<vx<<","<<vy<<","<<vz
	//	<<"Position uxyz: "<<ux<<","<<uy<<","<<uz
	//	<<G4endl;

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SetOptPhotonPolar()
{
	//G4double angle = G4UniformRand() * 360.0 * deg;
	G4double angle = 0.0 * deg;
	SetOptPhotonPolar(angle);
}

void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
	if(fParticleGun->GetParticleDefinition() !=
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
	G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
	G4ThreeVector product = normal.cross(kphoton);
	G4double modul2       = product * product;
	G4ThreeVector e_perpend(0., 0., 1.);
	
	if(modul2 > 0.)
	e_perpend = (1. / std::sqrt(modul2)) * product;
	G4ThreeVector e_paralle = e_perpend.cross(kphoton);
	G4ThreeVector polar =
	std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
	fParticleGun->SetParticlePolarization(polar);
}

void PrimaryGeneratorAction::SetRandomDirection(G4bool val)
{
	fRandomDirection = val;
}

void PrimaryGeneratorAction::SetRandomInLMO1(G4bool val)
{
	fRandomInLMO1 = val;
}










char readExtFile(const char* filename){
	FILE *p = fopen(filename,"r");
		if(p == NULL){
		return '0';
	}
	return fgetc(p);
}
