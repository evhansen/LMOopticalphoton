#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
//#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  //virtual void ConstructParticle();
  //virtual void ConstructProcess();

  virtual void SetCuts();

private:
  G4int ver;
};


#endif

