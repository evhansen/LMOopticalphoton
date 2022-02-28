
#ifndef LightDetectorSD_h
#define LightDetectorSD_h 1

#include "G4VSensitiveDetector.hh"

#include "LightDetectorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class LightDetectorSD : public G4VSensitiveDetector
{
  public:
    LightDetectorSD(const G4String& name,
                     const G4String& hitsCollectionName,
                     G4int nofCells);
    virtual ~LightDetectorSD();

    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    LightDetectorHitsCollection* fHitsCollection;
    G4int  fNofCells;
};

#endif
