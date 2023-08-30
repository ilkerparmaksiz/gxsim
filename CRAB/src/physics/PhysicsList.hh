//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4SystemOfUnits.hh"

#include "G4GenericIon.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
// particles
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "NESTProc.hh"


class G4VPhysicsConstructor;
class PhysicsListMessenger;
class DetectorConstruction;
class G4FastSimulationPhysics;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList : public G4VModularPhysicsList {
 public:
  PhysicsList();
  virtual ~PhysicsList();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);


  void SetLowerProductionLimit(G4double e){lowE=e;};

  void InitializePhysicsList(const G4String& name);
  void AddParametrisation();
  //   void ConstructParticle();
  //  void ConstructProcess() ; // uncommenting this makes my simulation suddenly have no e- physics, infinite track lengths
  
 private:
  void AddIonGasModels();
 

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double lowE;

  PhysicsListMessenger* pMessenger;
  G4FastSimulationPhysics* fastSimulationPhysics;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
