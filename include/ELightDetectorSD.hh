
#ifndef ELightDetectorSD_h
#define ELightDetectorSD_h 1

#include "G4VSensitiveDetector.hh"

#include "LightDetectorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class ELightDetectorSD : public G4VSensitiveDetector
{
  public:
    ELightDetectorSD(const G4String& name,
                     const G4String& hitsCollectionName,
                     G4int nofCells);
    virtual ~ELightDetectorSD();

    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    LightDetectorHitsCollection* fHitsCollection;
    G4int  fNofCells;
};

#endif
