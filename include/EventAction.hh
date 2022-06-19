#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "LightDetectorHit.hh"

#include "globals.hh"

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

private:
  LightDetectorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;

  void PrintEventStatistics(G4double absoEdep//, G4double absoTrackLength
                            // ,G4double gapEdep, G4double gapTrackLength
                          ) const;

  // data members
  G4int  fAbsHCID;
  // G4int  fGapHCID;
};

#endif
