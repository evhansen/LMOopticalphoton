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
fPhysVol(nullptr),
fPhysVolNum()
{}

LightDetectorHit::LightDetectorHit(G4VPhysicalVolume* pVol)
: fPhysVol(pVol)
{}

LightDetectorHit::LightDetectorHit(const LightDetectorHit& right)
: G4VHit()
{
	fEdep        = right.fEdep;
	fPhysVol     = right.fPhysVol;
	fPhysVolNum  = right.fPhysVolNum;
	fPos 	     = right.fPos;
	//iPos	     = right.iPos; 
}


LightDetectorHit::~LightDetectorHit() {}


/*
const LightDetectorHit& LightDetectorHit::operator=(const LightDetectorHit& right)
{
	fEdep        = right.fEdep;
	fPhysVol     = right.fPhysVol;
	fPhysVolNum  = right.fPhysVolNum;
	return *this;
}

G4bool LightDetectorHit::operator==(const LightDetectorHit& right) const
{
	return ( this == &right ) ? true : false;
}
*/

void LightDetectorHit::Print()
{
	G4cout
	<< "\nEdep: "
	<< std::setw(7) << G4BestUnit(fEdep,"Energy")
	<< "\nPhysVol: "
	<< std::setw(7) << fPhysVol
	<< "\nPhysVolNum: "
	<< std::setw(7) << fPhysVolNum
	<< "\nfPosx,y,z: "
	<< std::setw(7) << fPos.getX() << "," << fPos.getY() << "," << fPos.getZ()
	<< G4endl;
}
