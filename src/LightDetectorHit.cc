#include "LightDetectorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

G4ThreadLocal G4Allocator<LightDetectorHit>* LightDetectorHitAllocator = 0;

LightDetectorHit::LightDetectorHit()
: G4VHit(),
fEdep(0.),
fPhysVolNum(-1)
{}



LightDetectorHit::~LightDetectorHit() {}



void LightDetectorHit::Print()
{
	G4cout
	<< "\nEdep: "
	<< std::setw(7) << G4BestUnit(fEdep,"Energy")
	<< "\nPhysVolNum: "
	<< std::setw(7) << fPhysVolNum
	<< "\nfPosx,y,z: "
	<< std::setw(7) << fPos.getX() << "," << fPos.getY() << "," << fPos.getZ()
	<< G4endl;
}
