#include "SteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "../geometry/CRAB_CSG.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4StepPoint.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "GasBoxSD.hh"
#include "S2Photon.hh"
#include "config.h"
#ifdef With_Opticks
#include "G4CXOpticks.hh"
#include "SEvt.hh"
#include "U4.hh"
#include "G4Scintillation.hh"
#endif
#include "G4Alpha.hh"
namespace {G4Mutex Mutex= G4MUTEX_INITIALIZER;}
SteppingAction::SteppingAction(EventAction *eva) : fEventAction(eva),  ev_shift(0) {

  msg_ = new G4GenericMessenger(this, "/Action/SteppingAction/",
    "Control commands of the stepping action.");

  msg_->DeclareProperty("event_shift", ev_shift, "Set the event ID Shift number");
  
}



void SteppingAction::UserSteppingAction(const G4Step *aStep)

{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  const G4ThreeVector pos(aStep->GetPreStepPoint()->GetPosition());
  const G4ThreeVector tpos(aStep->GetPostStepPoint()->GetPosition());
  const G4ThreeVector sdir(aStep->GetPreStepPoint()->GetMomentumDirection());
  const G4ThreeVector tdir(aStep->GetPostStepPoint()->GetMomentumDirection());

  
  std::string startp("null");
  std::string endp("null");

  const G4StepPoint* endPoint = aStep->GetPostStepPoint();

  const G4VProcess* tprocess   = endPoint->GetProcessDefinedStep();
  const G4VProcess* sprocess   = aStep->GetPreStepPoint()->GetProcessDefinedStep();

  tprocess = aStep->GetPostStepPoint()->GetProcessDefinedStep();

  G4Track* track = aStep->GetTrack();
  G4double tID = track->GetTrackID();

  /*if (trackID != track->GetTrackID()){
    trackID =  track->GetTrackID();
    reflected = false;
  }*/

  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4int pID       = particle->GetPDGEncoding();
  G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();

  G4OpBoundaryProcess* boundary = 0;

  // if (!boundary &&  particle->GetParticleName() == "S2Photon") {
  if (!boundary){
      
      G4ProcessVector* pv = particle->GetProcessManager()->GetProcessList();
      for (size_t i=0; i<pv->size(); i++) {
          
          if ((*pv)[i]->GetProcessName() == "OpBoundary") {
              boundary = (G4OpBoundaryProcess*) (*pv)[i];
              
              /*// This is a reflection on the SS
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "GXe" && aStep->GetPostStepPoint()->GetMaterial()->GetName() == "Steel" && boundary->GetStatus() == SpikeReflection){
                  reflected = true;
               }

              // This gives the material that we just reflected off
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "Steel" && aStep->GetPostStepPoint()->GetMaterial()->GetName() == "GXe" && boundary->GetStatus() == StepTooSmall){
                Material_Store = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }

              // This gives the material that we just reflected off
              if(aStep->GetPostStepPoint()->GetMaterial()->GetName() == "MgF2" && boundary->GetStatus() == FresnelReflection){
                reflected = true;
                Material_Store = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }

              // Total internal reflection from MgF2
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "MgF2" && boundary->GetStatus() == TotalInternalReflection){
                reflected = true;
                Material_Store = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }*/
              if(boundary->GetStatus() == Detection  )  {
                  Material_Store = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
                  G4String detector_name = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

                  //std::cout<<"saving to file Camera" << std::endl;
                  analysisManager->FillNtupleDColumn(0,0, event+ev_shift);
                  analysisManager->FillNtupleDColumn(0,1, aStep->GetTrack()->GetParentID());
                  analysisManager->FillNtupleDColumn(0,2, time/ns);
                  analysisManager->FillNtupleDColumn(0,3, tpos[0]/mm); // Get Post step position
                  analysisManager->FillNtupleDColumn(0,4, tpos[1]/mm);
                  analysisManager->FillNtupleDColumn(0,5, tpos[2]/mm);
                  analysisManager->FillNtupleSColumn(0,6, detector_name.c_str());
                  analysisManager->AddNtupleRow(0);
              }
              break;
          }
      }
  }

    #ifdef With_Opticks

    if(aStep->GetTrack()->GetDefinition()==G4OpticalPhoton::Definition() and aStep->GetTrack()->GetDefinition()==S2Photon::Definition() ) {
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        return;
    }
    //if(particle==G4Alpha::Definition())  SEvt::AddTorchGenstep();
    G4SteppingManager * sMg=G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
    G4StepStatus stepStatus=sMg->GetfStepStatus();
    if(stepStatus!=fAtRestDoItProc){
        G4ProcessVector * PostStepProc=sMg->GetfPostStepDoItVector();
        size_t MaxSteps=sMg->GetMAXofPostStepLoops();
        for (int stp=0;stp<MaxSteps;stp++){
            if((*PostStepProc)[stp]->GetProcessName()=="Scintillation"){
                G4Scintillation *ScintProc= (G4Scintillation*) (*PostStepProc)[stp];
                G4int num_photons=ScintProc->GetNumPhotons();
                std::cout << "Scintilation "<< num_photons <<std::endl;

                if(num_photons>0){

                    G4MaterialPropertiesTable * MPT = track->GetMaterial()->GetMaterialPropertiesTable();
                    G4double t1,t2=0;
                    G4int singlets,triplets=0;
                    t1=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
                    t2=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
                    singlets= floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1)*num_photons);
                    triplets= ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2)*num_photons);
                    //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
                   if(singlets>0)
                        U4::CollectGenstep_DsG4Scintillation_r4695(track,aStep,singlets,0,t1);
                    if(triplets>0)
                        U4::CollectGenstep_DsG4Scintillation_r4695(track,aStep,triplets,1,t2);

                }

            }
        }
    }
    #endif



  G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
                             ->GetVolume()->GetLogicalVolume();

   
  G4TouchableHandle touch = endPoint->GetTouchableHandle();

  if (pID==11 && track->GetKineticEnergy()/keV>0.100 && (lVolume->GetName().find("GAS_")!=std::string::npos)) // don't count the thermale's, just G4 e's
    fEventAction->EDepPrim(aStep->GetTotalEnergyDeposit());


  /*// Here we select if we got an S1 photon or S2 photon
  // S1  = 1, S2 = 2
  G4int PhotonType = 2;
  if(track->GetParticleDefinition()==S2Photon::OpticalPhoton()){
    PhotonType = 2;
  }
  else if(track->GetParticleDefinition()==G4OpticalPhoton::OpticalPhoton()){
    PhotonType = 1;
  }
    */


}



  
