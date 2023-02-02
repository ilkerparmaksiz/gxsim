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

  if (pID==11 && track->GetKineticEnergy()/keV>0.100 && (lVolume->GetName().find("GAS")!=std::string::npos)) // don't count the thermale's, just G4 e's
    fEventAction->EDepPrim(aStep->GetTotalEnergyDeposit());
 
 

  
  if (eVolume)
    {
      eVname = eVolume->GetName();

      if (lVolume->GetName().find("camLogical")!=std::string::npos) // PMT
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
	  //	  analysisManager->FillNtupleSColumn(id,6, startp);

	  analysisManager->AddNtupleRow(id);

	}
      id = 4;

      if (lVolume->GetName().find("Lens")!=std::string::npos) //      
	//      if (/*eVname=="detectorLogical" && */ lVolume->GetName()=="pmtPhysical")
      	{
	  //	  std::cout << "SteppingAction: Stepping from  " << lVolume->GetName() <<  " into " << eVname << " Killing OpticalPhoton." << std::endl;
	  //std::cout << "SteppingAction: startp and endproc are " << startp << " and " << endp << std::endl;

	  analysisManager->FillNtupleDColumn(id,0, event);
	  analysisManager->FillNtupleDColumn(id,1, pID);
	  analysisManager->FillNtupleDColumn(id,2, time/ns);
	  analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
	  analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
	  analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);
	  //	  analysisManager->FillNtupleSColumn(id,6, startp);

	  /*
	  std::cout << "current, exiting Volumes" << lVolume->GetName() << ", " << eVolume->GetName() << std::endl;
	  std::cout << "Incoming lens dir:" << sdir[0]<<"," <<sdir[1] << "," <<sdir[2] <<" exiting dir " << eVname << " at: "<< tdir[0]  <<"," <<tdir[1] << "," <<tdir[2] << std::endl;
	  std::cout << "In/out lens pos:" << pos[0]<<"," <<pos[1] << "," <<pos[2] <<" exiting pos " << eVname << " at: "<< tpos[0]  <<"," <<tpos[1] << "," <<tpos[2] << std::endl;
	  std::cout << "In lens energy; trackID, startp, endp, trackstatus : " << tID << ", " << aStep->GetPreStepPoint()->GetTotalEnergy() << "; " << startp << ", " << endp << ", " << track->GetTrackStatus() <<std::endl;
	  */
	  
	  analysisManager->AddNtupleRow(id);
  
	  

	}


    

    }
  
  // if particle == thermale, opticalphoton and parent == primary and stepID==1, or trackID<=2
  // count the NEST e-s/photons into a class variable from the primary particle. Retrieve at EndEvent().
}



  
