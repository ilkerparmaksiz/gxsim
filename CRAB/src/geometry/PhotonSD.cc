//
// Created by argon on 2/6/24.
//

#include "PhotonSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
// project headers
#include "PhotonSD.hh"
#include "scuda.h"
#include "SEvt.hh"
#include "G4CXOpticks.hh"
#include "NP.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
PhotonSD::PhotonSD(G4String name)
        : G4VSensitiveDetector(name)
{
    G4String HCname = name + "_HC";
    collectionName.insert(HCname);
    G4cout << collectionName.size() << "   PhotonSD name:  " << name << " collection Name: " << HCname
           << G4endl;
    fHCID = -1;
    //  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
}

void PhotonSD::Initialize(G4HCofThisEvent* hce)
{
    fPhotonHitsCollection = new PhotonHitsCollection(SensitiveDetectorName, collectionName[0]);
    if(fHCID < 0)
    {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    hce->AddHitsCollection(fHCID, fPhotonHitsCollection);
    fOpticksHits=0;
    fGeant4Hits=0;
}

G4bool PhotonSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
    G4TouchableHistory* theTouchable =
            (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physVol = theTouchable->GetVolume();
    G4Track* theTrack          = aStep->GetTrack();
    // we only deal with optical Photons:



    if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
    {
        return false;
    }
    auto ana = G4AnalysisManager::Instance();
    auto run =G4RunManager::GetRunManager();

    G4double theEdep              = theTrack->GetTotalEnergy() / CLHEP::eV;
    const G4VProcess* thisProcess = theTrack->GetCreatorProcess();
    G4String processname;
    if(thisProcess != NULL)
        processname = thisProcess->GetProcessName();
    else
        processname = "No Process";
    int theCreationProcessid;
    if(processname == "Cerenkov")
    {
        theCreationProcessid = 0;
    }
    else if(processname == "Scintillation")
    {
        theCreationProcessid = 1;
    }
    else
    {
        theCreationProcessid = -1;
    }
    //std::cout << processname << " " << aStep->GetPostStepPoint()->GetPosition() <<std::endl;
    PhotonHit* newHit = new PhotonHit(0,theCreationProcessid, etolambda(theEdep), theTrack->GetGlobalTime(),
            aStep->GetPostStepPoint()->GetPosition(), aStep->GetPostStepPoint()->GetMomentumDirection(),
            aStep->GetPostStepPoint()->GetPolarization());
    fPhotonHitsCollection->insert(newHit);

    ana->FillNtupleIColumn(0,0,run->GetCurrentEvent()->GetEventID());
    ana->FillNtupleDColumn(0,1,theTrack->GetPosition()[0]/CLHEP::mm);
    ana->FillNtupleDColumn(0,2,theTrack->GetPosition()[1]/CLHEP::mm);
    ana->FillNtupleDColumn(0,3,theTrack->GetPosition()[2]/CLHEP::mm);
    ana->FillNtupleDColumn(0,4,theTrack->GetProperTime()/CLHEP::ns);
    ana->FillNtupleSColumn(0,5,physVol->GetName());
    ana->FillNtupleSColumn(0,6,processname);
    ana->AddNtupleRow(0);
    fGeant4Hits++;
     theTrack->SetTrackStatus(fStopAndKill);
    return true;
}

void PhotonSD::EndOfEvent(G4HCofThisEvent*)
{
    // if(ConfigurationManager::getInstance()->isEnable_verbose())
    //{
    G4int NbHits = fPhotonHitsCollection->entries();
    G4cout << " PhotonSD::EndOfEvent Number of PhotonHits:  " << NbHits << G4endl;
    //}
}
void PhotonSD::OpticksHits()
{
  SEvt* sev             = SEvt::Get_EGPU();
  unsigned int num_hits = sev->GetNumHit(0);
  auto ana=G4AnalysisManager::Instance();
  auto run= G4RunManager::GetRunManager();
  fOpticksHits=sev->GetNumHit(0);
  for(int idx = 0; idx < int(num_hits); idx++)
  {
    sphoton hit;
    sev->getHit(hit, idx);
    G4ThreeVector position     = G4ThreeVector(hit.pos.x, hit.pos.y, hit.pos.z);
    G4ThreeVector direction    = G4ThreeVector(hit.mom.x, hit.mom.y, hit.mom.z);
    G4ThreeVector polarization = G4ThreeVector(hit.pol.x, hit.pol.y, hit.pol.z);
    int theCreationProcessid;
    if(OpticksPhoton::HasCerenkovFlag(hit.flagmask))
    {
      theCreationProcessid = 0;
    }
    else if(OpticksPhoton::HasScintillationFlag(hit.flagmask))
    {
      theCreationProcessid = 1;
    }
    else
    {
      theCreationProcessid = -1;
    }
/*
    PhotonHit* newHit = new PhotonHit(hit.iindex, theCreationProcessid, hit.wavelength, hit.time,
                                      position, direction, polarization);
    fPhotonHitsCollection->insert(newHit);
    //

    G4cout << " Process ID: " << theCreationProcessid << " PhotonSD  pos.:" << hit.pos.x << "  "
           << hit.pos.y << "  "
           << "  " << hit.pos.z << "  mom.:  " << hit.mom.x << "  " << hit.mom.y << "  "
           << hit.mom.z << "  pol.:  " << hit.pol.x << "  "
           << "  " << hit.pol.y << "  " << hit.pol.z << " iiindex: " << hit.iindex << "  "
           << "  wavel.:  " << hit.wavelength << "  time:  " << hit.time
           << "  boundary flag:  " << hit.boundary_flag << "  identy:  " << hit.identity
           << "  orient_idx: " << hit.orient_idx << "  flagmask:  " << hit.flagmask << G4endl;

*/
      // Lets save the hits
      ana->FillNtupleIColumn(2,0,run->GetCurrentEvent()->GetEventID());
      ana->FillNtupleIColumn(2,1,hit.iindex);
      ana->FillNtupleFColumn(2,2,hit.pos.x);
      ana->FillNtupleFColumn(2,3,hit.pos.y);
      ana->FillNtupleFColumn(2,4,hit.pos.z);
      ana->FillNtupleFColumn(2,5,hit.time);
      ana->FillNtupleFColumn(2,6,hit.mom.x);
      ana->FillNtupleFColumn(2,7,hit.mom.y);
      ana->FillNtupleFColumn(2,8,hit.mom.z);
      ana->FillNtupleFColumn(2,9,hit.pol.x);
      ana->FillNtupleFColumn(2,10,hit.pol.y);
      ana->FillNtupleFColumn(2,11,hit.pol.z);
      ana->FillNtupleFColumn(2,12,hit.wavelength);
      ana->AddNtupleRow(2);
  }
}

