#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "S2Photon.hh"
//#include "G4QStoppingPhysics.hh"//not in g4 10
#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
//#include "G4HadronQElasticPhysics.hh"	//not in g4 10
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ProcessManager.hh"

#include "DegradModel.hh"
#include "GarfieldVUVPhotonModel.hh"



#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4OpticalPhysics.hh"

/// Below is not found in my 10.7.2 build. Try commenting out to see if we get away w it. EC, 12-Oct-2021.
/////#include "G4OpticalProcessIndex.hh"

#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"

#include "PhysListEmStandard.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4PAIPhotModel.hh"
#include "G4PAIModel.hh"

#include "G4StepLimiter.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"

#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "NESTS1Photon.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4GlobalFastSimulationManager.hh"

#ifdef theParticleIterator
#undef theParticleIterator
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
    : G4VModularPhysicsList(), lowE(-1){
  G4LossTableManager::Instance();
  defaultCutValue = 10. * um;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForPositron = defaultCutValue;
  // Particles
  S2Photon::Definition();
  NEST::NESTThermalElectron::Definition();
  NESTS1Photon::Definition();
  pMessenger = new PhysicsListMessenger(this);


  SetVerboseLevel(0);

  RegisterPhysics(new G4EmLivermorePhysics(1));

  // Add General Decay
  RegisterPhysics(new G4DecayPhysics());

  // Add final state interactions following radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  RegisterPhysics(new G4StepLimiterPhysics());
  
  //This is needed to notify Geant4 that the G4FastSimulationModel is to be used as a possible physics process
  fastSimulationPhysics = new G4FastSimulationPhysics("fastSimPhys");
  RegisterPhysics(fastSimulationPhysics);

  
  //  RegisterPhysics(new G4OpticalPhysics());

  
  std::cout << "Finished PhysicsList constructor" << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {
  delete pMessenger;
  G4cout << "Deleting PhysicsList" << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::InitializePhysicsList(const G4String& name) {
  if (verboseLevel > 1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == "local") {
    ReplacePhysics(new PhysListEmStandard(name));
  } else if (name == "emstandard_opt0") {
    ReplacePhysics(new G4EmStandardPhysics(1));
  } else if (name == "emstandard_opt1") {
    ReplacePhysics(new G4EmStandardPhysics_option1());
  } else if (name == "emstandard_opt2") {
    ReplacePhysics(new G4EmStandardPhysics_option2());
  } else if (name == "emstandard_opt3") {
    ReplacePhysics(new G4EmStandardPhysics_option3());
  } else if (name == "emlivermore") {
    ReplacePhysics(new G4EmLivermorePhysics());
    std::cout << "PhysicsList::InitializePhysics initialised EMLivermore." << std::endl;
  } else if (name == "empenelope") {
    ReplacePhysics(new G4EmPenelopePhysics()); 
  } else if (name == "ionGasModels") {
    ReplacePhysics(new G4EmStandardPhysics(1));
    AddIonGasModels();
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts() {
  if (verboseLevel > 0) {
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
  }
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  G4cout << "The lower energy production cut is set to: " << lowE << G4endl;
  if(lowE>0)
      G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowE, 100. * MeV);
  
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e+"));
  
  if (region) {
    region->SetProductionCuts(cuts);
  }

  G4cout << "The lower energy production cut after call to GetProductionCutsTable()->SetEnergyRange(lowE, 100. * MeV) is now set to: " << lowE << G4endl;
  if (verboseLevel > 0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut) {
  cutForGamma = cut;
  G4cout << "Setting cut for gamma " << G4BestUnit(cut, "Length") << G4endl;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut) {
  cutForElectron = cut;
  G4cout << "Setting cut for electron " << G4BestUnit(cut, "Length") << G4endl;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut) {
  cutForPositron = cut;
  G4cout << "Setting cut for positron " << G4BestUnit(cut, "Length") << G4endl;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddIonGasModels() {
  G4EmConfigurator* em_config =
      G4LossTableManager::Instance()->EmConfigurator();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String partname = particle->GetParticleName();
    if (partname == "alpha" || partname == "He3" || partname == "GenericIon") {
      G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
      G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
      G4double eth = 2. * MeV * particle->GetPDGMass() / proton_mass_c2;
      em_config->SetExtraEmModel(partname, "ionIoni", mod1, "", 0.0, eth,
                                 new G4IonFluctuations());
      em_config->SetExtraEmModel(partname, "ionIoni", mod2, "", eth, 100 * TeV,
                                 new G4UniversalFluctuation());
    }
  }
}

//This activates the G4FastSimulationPhysics for all particles and should be called by the user in the macro before '/run/initialize' (command: '/MWGPC/phys/AddParametrisation')
void PhysicsList::AddParametrisation() {   
    theParticleTable->GetIterator()->reset();
    while ((*theParticleTable->GetIterator())()) {
        G4String particleName = theParticleTable->GetIterator()->value()->GetParticleName();
        fastSimulationPhysics->ActivateFastSimulation(particleName);
    }
}




/*
void PhysicsList::ConstructProcess()
{
  


    
  gasNESTdet* gndet = new gasNESTdet();
  NEST::NESTcalc* calcNEST = new NEST::NESTcalc(gndet);

  //  theNEST2ScintillationProcess->SetDetailedSecondaries(true);
  // theNEST2ScintillationProcess->SetStackElectrons(true);
OA  

  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic, calcNEST, gndet); //gndet);
    
    if (theNEST2ScintillationProcess->IsApplicable(*particle) ) {
      std::cout << "PhysicsList::ConstructProcess(): particleName, pmanager  " << particleName << ", " << pmanager << "." << std::endl;
      std::cout << "ordDefault, ordInActive " << ordDefault << ", " << ordInActive  << std::endl;
      pmanager->AddProcess(theNEST2ScintillationProcess, ordDefault + 1, ordInActive, ordDefault + 1);
    }


    if ( !(particleName.find("e-")!=std::string::npos  || particleName.find("opticalphoton")!=std::string::npos ) )
      continue;

    std::cout << "PhysicsList::ConstructProcess(): e-/optphot in particleIterator" << std::endl;
    //    G4Scintillation* fScintillationProcess = new G4Scintillation("Scintillation");
    //    G4OpBoundaryProcess* fBoundaryProcess = new G4OpBoundaryProcess();

    if (particleName == "opticalphoton" && pmanager) {
      //      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      //     pmanager->AddDiscreteProcess(fAbsorptionProcess);
      //      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
      //      pmanager->AddDiscreteProcess(fTheWLSProcess);
      //      pmanager->AddDiscreteProcess(fBoundaryProcess);
    }
    
    
  }



    
}
*/



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
