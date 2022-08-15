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
#include "G4ParticleDefinition.hh"

/// LightDetector hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class LightDetectorHit : public G4VHit
{
  public:
    LightDetectorHit();
    
    virtual ~LightDetectorHit();

    // operators
    const LightDetectorHit& operator=(const LightDetectorHit& rLDH);
    /*{
    	
	fEdep = rLDH.fEdep;
	fPhysVolNum = rLDH.fPhysVolNum;
    	fPos = rLDH.fPos;
	//iPos = rLDH.iPos;
	
	return *this;
    }*/
    	
    G4bool operator==(const LightDetectorHit& rLDH) const {
	return (this==&rLDH) ? true : false;
    }

    inline void* operator new(size_t);
    inline void  operator delete(void*);


    // methods from base class
    virtual void Draw() {}
    virtual void Print();


    // get and set methods
    G4double GetEdep() const {return fEdep;}
    inline void SetEdep(G4double de) { fEdep = de; }
    inline void AddEdep(G4double de) { fEdep += de; }

    G4double GetPhysVolNum() const {return fPhysVolNum;}
    inline void SetPhysVolNum(G4int plv) { fPhysVolNum = plv; }

    G4double GetPDG() const {return PDG;}
    inline void SetPDG(G4int plv) { PDG = plv; }

    const G4ParticleDefinition* GetPDef() const {return fPDef;}
    inline void SetPDef(const G4ParticleDefinition* pid) { fPDef = pid; }

    G4double GetEnergy() const {return Energy;}
    inline void SetEnergy(G4double eng) { Energy = eng; }

    G4double GetxfPos() const {return fPos.getX();}
    G4double GetyfPos() const {return fPos.getY();}
    G4double GetzfPos() const {return fPos.getZ();}
    inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }

    G4double GetxiPos() const {return iPos.getX();}
    G4double GetyiPos() const {return iPos.getY();}
    G4double GetziPos() const {return iPos.getZ();}
    inline void SetiPos(G4ThreeVector xyz) { iPos = xyz; }

    G4double GetxiMom() const {return iMom.getX();}
    G4double GetyiMom() const {return iMom.getY();}
    G4double GetziMom() const {return iMom.getZ();}
    inline void SetiMom(G4ThreeVector xyz) { iMom = xyz; }
 
    G4double GetxfMom() const {return fMom.getX();}
    G4double GetyfMom() const {return fMom.getY();}
    G4double GetzfMom() const {return fMom.getZ();}
    inline void SetfMom(G4ThreeVector xyz) { fMom = xyz; }
    

    G4int GetEventID() const {return fEventID;}
    G4int GetTrackID() const {return fTrackID;}
    G4int GetParentID() const {return fParentID;}
    inline void SetTrackID(G4int tid) { fTrackID = tid; }
    inline void SetParentID(G4int pid) { fParentID = pid; }
    inline void SetEventID(G4int eid) { fEventID = eid; }

  private:
    
    G4double Energy;
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    
    G4int fPhysVolNum;
    G4int PDG;

    G4ThreeVector fPos;
    G4ThreeVector iPos;

    G4ThreeVector fMom;
    G4ThreeVector iMom;

    G4int fEventID;
    G4int fTrackID;
    G4int fParentID; 
    const G4ParticleDefinition* fPDef;
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


#endif
