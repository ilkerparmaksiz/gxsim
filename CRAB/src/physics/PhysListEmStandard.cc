//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: PhysListEmStandard.cc,v 1.24 2009-11-15 22:10:03 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"
#include "G4SystemOfUnits.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4KleinNishinaModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4EmParameters.hh"
#include "G4MscStepLimitType.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4SystemOfUnits.hh"

#include "G4EmModelActivator.hh"

#include "G4FastSimulationManagerProcess.hh"

#include "OpWLS.hh"
#include "OpMieHG.hh"
#include "OpRayleigh.hh"
#include "OpBoundaryProcess.hh"
#include "OpWLS2.hh"

#include "gasNESTdet.hh"
#include "G4/NESTProc.hh"
#include "OpAbsorption.hh"

#include "S2Photon.hh"
#ifdef theParticleIterator
#undef theParticleIterator
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
: G4VPhysicsConstructor(name){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard() {G4cout << "Deleting PhysListEmStandard" << G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess() {
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    
    // Add standard EM Processes
    //
    G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
    theParticleIterator->reset();
    while ((*theParticleIterator)() ) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String particleName = particle->GetParticleName();
        
        if (particleName == "gamma") {
            ////ph->RegisterProcess(new G4RayleighScattering, particle);
            ph->RegisterProcess(new G4PhotoElectricEffect, particle);
            G4ComptonScattering* cs = new G4ComptonScattering;
            cs->SetEmModel(new G4KleinNishinaModel());
            ph->RegisterProcess(cs, particle);
            ph->RegisterProcess(new G4GammaConversion, particle);
            
        } else if (particleName == "e-") {
            ph->RegisterProcess(new G4eMultipleScattering(), particle);
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.1, 100 * um);
            eIoni->SetDEDXBinning(12 * 10);
            ph->RegisterProcess(eIoni, particle);
            ph->RegisterProcess(new G4eBremsstrahlung(), particle);
            
        } else if (particleName == "e+") {
            ph->RegisterProcess(new G4eMultipleScattering(), particle);
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.1, 100 * um);
            eIoni->SetDEDXBinning(12 * 10);
            ph->RegisterProcess(eIoni, particle);
            ph->RegisterProcess(new G4eBremsstrahlung(), particle);
            ph->RegisterProcess(new G4eplusAnnihilation(), particle);
            
        } else if (particleName == "mu+" || particleName == "mu-") {
            ph->RegisterProcess(new G4MuMultipleScattering(), particle);
            G4MuIonisation* muIoni = new G4MuIonisation();
            muIoni->SetStepFunction(0.1, 50 * um);
            ph->RegisterProcess(muIoni, particle);
            ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
            ph->RegisterProcess(new G4MuPairProduction(), particle);
            
        } else if (particleName == "proton" || particleName == "pi-" ||
                   particleName == "pi+") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            G4hIonisation* hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.1, 20 * um);
            ph->RegisterProcess(hIoni, particle);
            ph->RegisterProcess(new G4hBremsstrahlung(), particle);
            ph->RegisterProcess(new G4hPairProduction(), particle);
            
        } else if (particleName == "alpha" || particleName == "He3") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetStepFunction(0.1, 1 * um);
            ionIoni->SetDEDXBinning(12 * 10);
            ph->RegisterProcess(ionIoni, particle);
            ph->RegisterProcess(new G4NuclearStopping(), particle);
            
        } else if (particleName == "GenericIon") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            G4ionIonisation* ionIoni = new G4ionIonisation();
            ionIoni->SetEmModel(new G4IonParametrisedLossModel());
            ionIoni->SetStepFunction(0.1, 1 * um);
            ph->RegisterProcess(ionIoni, particle);
            ph->RegisterProcess(new G4NuclearStopping(), particle);
            
        } else if ((!particle->IsShortLived()) &&
                   (particle->GetPDGCharge() != 0.0) &&
                   (particle->GetParticleName() != "chargedgeantino")) {
            // all others charged particles except geantino
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            ph->RegisterProcess(new G4hIonisation(), particle);
        }/*else if(particleName=="S2Photon"){
            ph->RegisterProcess(new OpBoundaryProcess,particle);
            ph->RegisterProcess(new OpAbsorption,particle);
        }*/
    }
    
    // Em options
    //
    // Main options and setting parameters are shown here.
    // Several of them have default values.
    //
    G4EmParameters *emOptions= G4EmParameters::Instance();
    
    // physics tables
    //
    emOptions->SetMinEnergy(10 * eV);      // default 100 eV
    emOptions->SetMaxEnergy(10 * TeV);     // default 100 TeV
    // emOptions.SetDEDXBinning(12 * 10);    // default=12*7
    // emOptions.SetLambdaBinning(12 * 10);  // default=12*7
    
    // multiple coulomb scattering
    //
    emOptions->SetMscStepLimitType(fUseSafety); // default  // default
    
    // Deexcitation
    //
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    de->SetFluo(true);
    de->SetAuger(false);
    de->SetPIXE(true);
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
    
    G4EmModelActivator mact(GetPhysicsName());


  auto particleIteratorP=GetParticleIterator();
  particleIteratorP->reset();
  std::cout  << "PhysicsListEMStandard::ConstructProcess() pit is "  << particleIteratorP << std::endl;

  gasNESTdet* gndet = new gasNESTdet();
  // std::shared_ptr<gasNESTdet> gndet(new gasNESTdet());
  NEST::NESTcalc* calcNEST = new NEST::NESTcalc(gndet);  

  NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic, calcNEST, gndet); //gndet);
  theNEST2ScintillationProcess->SetDetailedSecondaries(true);
  theNEST2ScintillationProcess->SetStackElectrons(true);

  while( (*particleIteratorP)() ){
    G4ParticleDefinition* particle = particleIteratorP->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // td::cout << "PhysicsListEMStandard::ConstructProcess(): pname, pmanager are " << particleName << ", " << pmanager << std::endl;
    // if ( !( particleName.find("e-")!=std::string::npos || particleName.find("alpha")!=std::string::npos  || particleName.find("opticalphoton")!=std::string::npos ) )
    //   continue;
      if (pmanager) {
        if (theNEST2ScintillationProcess->IsApplicable(*particle) && pmanager) {
        std::cout << "PhysicsList::InitialisePhysics(): particleName, pmanager  " << particleName << ", " << pmanager << "." << std::endl;
        std::cout << "ordDefault, ordInActive " << ordDefault << ", " << ordInActive  << std::endl;
        // This needs to be called because Nest is not aware of the S2Photon
        if(particleName!="S2Photon") pmanager->AddProcess(theNEST2ScintillationProcess, ordDefault + 1, ordInActive, ordDefault + 1);
      }
        OpBoundaryProcess* fBoundaryProcess = new OpBoundaryProcess();
        OpAbsorption* fAbsorptionProcess = new OpAbsorption();
        OpWLS* fTheWLSProcess = new OpWLS();
        OpMieHG* fOpMieHG = new OpMieHG();
        OpRayleigh* fTOpRayleigh = new OpRayleigh();
        OpWLS2* fOpWLS2 = new OpWLS2();

    if (((particleName == "opticalphoton") || particleName=="S2Photon") && pmanager) {
            G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
           pmanager->AddDiscreteProcess(fAbsorptionProcess);
           pmanager->AddDiscreteProcess(fBoundaryProcess);
	       //pmanager->AddDiscreteProcess(fTOpRayleigh);
           //pmanager->AddDiscreteProcess(fTheWLSProcess);
           //pmanager->AddDiscreteProcess(fOpMieHG);
           //pmanager->AddDiscreteProcess(fOpWLS2);
    }

      
      //      std::cout << "PhysicsList::InitialisePhysicsList(): e-/optphot in particleIteratorP" << std::endl;
      std::cout << "PhysicsListEMStandard()::ConstructProcess() process list (of length) " << pmanager->GetProcessList()->size() << " for " << particleName << " is: " << std::endl;
    //    for (const auto  proc : pmanager->GetProcessList())
      for (short int ii = 0; ii<(short int)(pmanager->GetProcessList()->size());++ii)
	{
	  std::cout << (*pmanager->GetProcessList())[ii]->GetProcessName() << std::endl;
	}
    }
  }


    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

