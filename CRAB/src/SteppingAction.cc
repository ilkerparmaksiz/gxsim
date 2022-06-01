#include "SteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "DetectorConstruction.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4StepPoint.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "GasBoxSD.hh"



SteppingAction::SteppingAction(){

}


void SteppingAction::UserSteppingAction(const G4Step *aStep) {

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  const G4ThreeVector pos(aStep->GetPreStepPoint()->GetPosition());
  const G4ThreeVector tpos(aStep->GetPostStepPoint()->GetPosition());

  std::string startp("null");
  std::string endp("null");

  const G4StepPoint* endPoint = aStep->GetPostStepPoint();

  const G4VProcess* tprocess   = endPoint->GetProcessDefinedStep();
  const G4VProcess* sprocess   = aStep->GetPreStepPoint()->GetProcessDefinedStep();

  tprocess = aStep->GetPostStepPoint()->GetProcessDefinedStep();

  G4Track* track = aStep->GetTrack();
  G4double tID = track->GetTrackID();

  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4int pID       = particle->GetPDGEncoding();
  G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();

  G4int id(0);

  if (sprocess)
      startp = sprocess->GetProcessName();
  if (tprocess)
      endp = tprocess->GetProcessName();

  
  G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
                             ->GetVolume()->GetLogicalVolume();

   
  G4TouchableHandle touch = endPoint->GetTouchableHandle();
  G4VPhysicalVolume* eVolume = touch->GetVolume();
  G4String eVname("null");
  if (eVolume)
    {
      eVname = eVolume->GetName();

      if (lVolume->GetName().find("pmt")!=std::string::npos)
	//      if (/*eVname=="detectorLogical" && */ lVolume->GetName()=="pmtPhysical")
      	{
	  //	  std::cout << "SteppingAction: Stepping from  " << lVolume->GetName() <<  " into " << eVname << " Killing OpticalPhoton." << std::endl;
	  //std::cout << "SteppingAction: startp and endproc are " << startp << " and " << endp << std::endl;
	  track->SetTrackStatus(fStopAndKill);
	  analysisManager->FillNtupleDColumn(id,0, event);
	  analysisManager->FillNtupleDColumn(id,1, pID);
	  analysisManager->FillNtupleDColumn(id,2, time/ns);
	  analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
	  analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
	  analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);
	  analysisManager->AddNtupleRow(id);

	}
    }

  id = 1;
  if (tpos[1]/cm <9.  && track->GetCurrentStepNumber()==1) // 9-> 2*det->HalfZ
    {
      //      std::cout << "Stepping bulk, num, z is " << track->GetCurrentStepNumber() << ", "<< tpos[1]/cm << std::endl;
	  analysisManager->FillNtupleDColumn(id,0, event);
	  analysisManager->FillNtupleDColumn(id,1, pID);
	  analysisManager->FillNtupleDColumn(id,2, time/ns);
	  analysisManager->FillNtupleDColumn(id,3, tpos[0]/mm);
	  analysisManager->FillNtupleDColumn(id,4, tpos[1]/mm);
	  analysisManager->FillNtupleDColumn(id,5, tpos[2]/mm);
	  analysisManager->AddNtupleRow(id);
    }

  id = 2;
  if (tpos[1]/cm >9. && track->GetCurrentStepNumber()==1) // 9-> 2*det->HalfZ
    {
      //          std::cout << "Stepping LEM, num, z is " << track->GetCurrentStepNumber() << ", "<< tpos[1]/cm << std::endl;
	  analysisManager->FillNtupleDColumn(id,0, event);
	  analysisManager->FillNtupleDColumn(id,1, pID);
	  analysisManager->FillNtupleDColumn(id,2, time/ns);
	  analysisManager->FillNtupleDColumn(id,3, tpos[0]/mm);
	  analysisManager->FillNtupleDColumn(id,4, tpos[1]/mm);
	  analysisManager->FillNtupleDColumn(id,5, tpos[2]/mm);
	  analysisManager->FillNtupleDColumn(id,6, track->GetCurrentStepNumber());
	  analysisManager->AddNtupleRow(id);
    }



  
  
  // if particle == thermale, opticalphoton and parent == primary and stepID==1, or trackID<=2
  // count the NEST e-s/photons into a class variable from the primary particle. Retrieve at EndEvent().
}
