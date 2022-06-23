/// \file LightDetectorHit.hh
/// \brief Definition of the LightDetectorHit class

#ifndef LightDetectorHit_h
#define LightDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4VPhysicalVolume.hh"

/// LightDetector hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class LightDetectorHit : public G4VHit
{
  public:
    LightDetectorHit();
    //LightDetectorHit(G4VPhysicalVolume* pVol);
    //LightDetectorHit(const LightDetectorHit& right);
    
    virtual ~LightDetectorHit();

    // operators
    const LightDetectorHit& operator=(const LightDetectorHit& rLDH){
    	
	fEdep = rLDH.fEdep;
	fPhysVolNum = rLDH.fPhysVolNum;
    	fPos = rLDH.fPos;
	//iPos = rLDH.iPos;
	
	return *this;
    }
    	
    G4bool operator==(const LightDetectorHit& rLDH) const {
	return (this==&rLDH) ? true : false;
    }

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de);

    // get methods
    G4double GetEdep() const;
    inline void SetEdep(G4double de) { fEdep = de; }
    inline void AddEdep(G4double de) { fEdep += de; }

    G4double GetPhysVolNum() const;
    inline void SetPhysVolNum(G4int plv) { fPhysVolNum = plv; }

    G4double GetxfPos() const;
    G4double GetyfPos() const;
    G4double GetzfPos() const;
    inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    

    G4int GetTrackID() const;
    G4int GetParentID() const;
    inline void SetTrackID(G4int tid) {fTrackID = tid;}
    inline void SetParentID(G4int pid) {fParentID = pid;}

  private:
    
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4int fPhysVolNum;
    G4ThreeVector fPos;
    G4int fTrackID;
    G4int fParentID; 
};


using LightDetectorHitsCollection = G4THitsCollection<LightDetectorHit>;

extern G4ThreadLocal G4Allocator<LightDetectorHit>* LightDetectorHitAllocator;


inline void* LightDetectorHit::operator new(size_t)
{
  if (!LightDetectorHitAllocator) {
    LightDetectorHitAllocator = new G4Allocator<LightDetectorHit>;
  }
  void *hit;
  hit = (void *) LightDetectorHitAllocator->MallocSingle();
  return hit;
}

inline void LightDetectorHit::operator delete(void *hit)
{
  if (!LightDetectorHitAllocator) {
    LightDetectorHitAllocator = new G4Allocator<LightDetectorHit>;
  }
  LightDetectorHitAllocator->FreeSingle((LightDetectorHit*) hit);
}

inline void LightDetectorHit::Add(G4double de) {
  fEdep += de;
}

inline G4double LightDetectorHit::GetEdep() const {
  return fEdep;
}

inline G4double LightDetectorHit::GetPhysVolNum() const {
  return fPhysVolNum;
}

inline G4double LightDetectorHit::GetxfPos() const {
	return fPos.getX();
}

inline G4double LightDetectorHit::GetyfPos() const {
	return fPos.getY();
}

inline G4double LightDetectorHit::GetzfPos() const {
	return fPos.getZ();
}


inline G4int LightDetectorHit::GetTrackID() const {
	return fTrackID;
}


inline G4int LightDetectorHit::GetParentID() const {
	return fParentID;
}

#endif
