//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file optical/OpNovice2/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
  : G4VUserPrimaryGeneratorAction() , fParticleGun(0), fDetector(det)
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);

  fRandomDirection = false;
  fPolarized       = false;
  fPolarization    = 0.;
  // default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e+");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0 * ns);
  fParticleGun->SetParticlePosition(
    G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fParticleGun->SetParticleEnergy(500.0 * keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fRandomDirection)
  {
	  G4double cosTheta = 2*G4UniformRand() - 1.;
	  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
	  G4double phi = CLHEP::twopi*G4UniformRand();
	  G4double vx = sinTheta*std::cos(phi);
	  G4double vy = sinTheta*std::sin(phi);
	  G4double vz = cosTheta;
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));

	/*G4double theta = CLHEP::halfpi * G4UniformRand();
	G4double phi   = CLHEP::twopi * G4UniformRand();
	G4double x     = std::cos(theta);
	G4double y     = std::sin(theta) * std::sin(phi);
	G4double z     = std::sin(theta) * std::cos(phi);
	G4ThreeVector dir(x, y, z);
	fParticleGun->SetParticleMomentumDirection(dir); */
  }
  if(fRandomInLMO1)
  {
	  G4double Rmax = fDetector->GetLMOXSize() ;
	  G4double SepDist = fDetector->GetLMOSep() ;
	  G4double R2 = Rmax*2.0;
	  G4double ux = G4UniformRand()*R2 - Rmax; //+ 0.5*SepDist ;
	  G4double uy = G4UniformRand()*R2 - Rmax; //+ 0.5*SepDist ;
	  G4double uz = G4UniformRand()*R2 - Rmax ;

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
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::SetRandomDirection(G4bool val)
{
  fRandomDirection = val;
}

void PrimaryGeneratorAction::SetRandomInLMO1(G4bool val)
{
  fRandomInLMO1 = val;
}






// /// \file PrimaryGeneratorAction.cc
// /// \brief Implementation of the PrimaryGeneratorAction class
//
// #include "PrimaryGeneratorAction.hh"
// #include "PrimaryGeneratorMessenger.hh"
// #include "DetectorConstruction.hh"
//
//
// #include "G4RunManager.hh"
// #include "G4LogicalVolumeStore.hh"
// #include "G4LogicalVolume.hh"
// #include "G4Box.hh"
// #include "G4Event.hh"
// #include "G4ParticleGun.hh"
// #include "G4ParticleTable.hh"
// #include "G4ParticleDefinition.hh"
// #include "G4SystemOfUnits.hh"
// #include "Randomize.hh"
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// PrimaryGeneratorAction::PrimaryGeneratorAction()
//  : G4VUserPrimaryGeneratorAction(),
//    fParticleGun(nullptr),
//    fDetector(det)
// {
//   G4int nofParticles = 1;
//   fParticleGun = new G4ParticleGun(nofParticles);
//   fGunMessenger = new PrimaryGeneratorMessenger(this);
//
//   fRandomDirection = false;
//   fPolarized       = false;
//   fPolarization    = 0.;
//
//   // default particle kinematic
//   //
//   auto particleDefinition
//     = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
//   fParticleGun->SetParticleDefinition(particleDefinition);
//   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
//   fParticleGun->SetParticleEnergy(3*eV);
// }
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// PrimaryGeneratorAction::~PrimaryGeneratorAction()
// {
//   delete fParticleGun;
//   delete fGunMessenger
// }
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
// {
//   // This function is called at the begining of event
//
//   // In order to avoid dependence of PrimaryGeneratorAction
//   // on DetectorConstruction class we get world volume
//   // from G4LogicalVolumeStore
//   //
//
//   if(fRandomDirection)
//   {
//     G4double cosTheta = 2*G4UniformRand() - 1.;
//     G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
//     G4double phi = CLHEP::twopi*G4UniformRand();
//     G4double vx = sinTheta*std::cos(phi);
//     G4double vy = sinTheta*std::sin(phi);
//     G4double vz = cosTheta;
//     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
//   }
//
//   if(fRandomInLMO1)
//   {
//     G4double Rmax = fDetector->GetLMOXSize() ;
//     G4double SepDist = fDetector->GetLMOSep() ;
//     G4double R2 = Rmax*2.0;
//     G4double ux = G4UniformRand()*R2 + 0.5*SepDist ;
//     G4double uy = G4UniformRand()*R2 + 0.5*SepDist ;
//     G4double uz = G4UniformRand()*R2 - Rmax ;
//
//     fParticleGun->SetParticlePosition(G4ThreeVector(ux,uy,uz));
//   }
//
//   if(fParticleGun->GetParticleDefinition() ==
//    G4OpticalPhoton::OpticalPhotonDefinition())
// {
//   if(fPolarized)
//     SetOptPhotonPolar(fPolarization);
//   else
//     SetOptPhotonPolar();
// }
// fParticleGun->GeneratePrimaryVertex(anEvent);
// }
//
// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
