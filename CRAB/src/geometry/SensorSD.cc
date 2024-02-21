// ----------------------------------------------------------------------------
// nexus | SensorSD.cc
//
// This class is the sensitive detector that allows for the registration
// of the charge detected by a photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "SensorSD.hh"
#include "S2Photon.hh"
#include <G4OpticalPhoton.hh>
#include <G4SDManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4RunManager.hh>
#include <G4AnalysisManager.hh>
#include "config.h"
#ifdef With_Opticks
#include "SEvt.hh"
#include "sphoton.h"
#endif
namespace sensorsd {


  SensorSD::SensorSD(G4String sdname):
    G4VSensitiveDetector(sdname),
    naming_order_(0), sensor_depth_(0), mother_depth_(0)
  {
    // Register the name of the collection of hits
    collectionName.insert(GetCollectionUniqueName());
  }



  SensorSD::~SensorSD()
  {
  }



  G4String SensorSD::GetCollectionUniqueName()
  {
    return "SensorHitsCollection";
  }



  void SensorSD::Initialize(G4HCofThisEvent* HCE)
  {
    // Create a new collection of PMT hits
    HC_ = new SensorHitsCollection(this->GetName(), this->GetCollectionName(0));

    G4int HCID = G4SDManager::GetSDMpointer()->
      GetCollectionID(this->GetName()+"/"+this->GetCollectionName(0));

    HCE->AddHitsCollection(HCID, HC_);
  }



  G4bool SensorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
  {
    // Check whether the track is an optical photon
    G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
    if (pdef != G4OpticalPhoton::Definition() || pdef!=S2Photon::Definition()) return false;

    //GEANT4Hits(step);
    std::cout <<"-------- Hits ---------" <<std::endl;
    const G4VTouchable* touchable =
      step->GetPostStepPoint()->GetTouchable();

    G4int pmt_id = FindPmtID(touchable);

    sensorhit::SensorHit* hit = 0;
    for (size_t i=0; i<HC_->entries(); i++) {
      if ((*HC_)[i]->GetPmtID() == pmt_id) {
        hit = (*HC_)[i];
        break;
      }
    }

    // If no hit associated to this sensor exists already,
    // create it and set main properties
    if (!hit) {
      hit = new sensorhit::SensorHit();
      hit->SetPmtID(pmt_id);
      hit->SetBinSize(timebinning_);
      hit->SetPosition(touchable->GetTranslation());
      HC_->insert(hit);
    }

    G4double time = step->GetPostStepPoint()->GetGlobalTime();
    hit->Fill(time);

    return true;
  }



  G4int SensorSD::FindPmtID(const G4VTouchable* touchable)
  {
    G4int pmtid = touchable->GetCopyNumber(sensor_depth_);
    if (naming_order_ != 0) {
      G4int motherid = touchable->GetCopyNumber(mother_depth_);
      pmtid = naming_order_ * motherid + pmtid;
    }
    return pmtid;
  }

  void SensorSD::OpticksHits()
    {
#ifdef With_Opticks
        SEvt* sev             = SEvt::Get_EGPU();

        unsigned int num_hits = sev->GetNumHit(0);

        auto ana=G4AnalysisManager::Instance();
        auto run= G4RunManager::GetRunManager();
        sev->GetNumHit(0);
        int id=8;

        std::cout<< "Name is " << this->GetName() <<std::endl;
        for(int idx = 0; idx < int(num_hits); idx++)
        {
            sphoton hit;


            //std::cout << hit.descDetail()  << " id  " <<id <<std::endl;

            sev->getHit(hit, idx);
            if(hit.boundary()==29) id=7;
            else id=8;
            G4ThreeVector position     = G4ThreeVector(hit.pos.x, hit.pos.y, hit.pos.z);
            G4ThreeVector direction    = G4ThreeVector(hit.mom.x, hit.mom.y, hit.mom.z);
            G4ThreeVector polarization = G4ThreeVector(hit.pol.x, hit.pol.y, hit.pol.z);


            // Lets save the hits
            ana->FillNtupleIColumn(id,0,run->GetCurrentEvent()->GetEventID());
            ana->FillNtupleIColumn(id,1,hit.iindex);
            ana->FillNtupleFColumn(id,2,hit.pos.x);
            ana->FillNtupleFColumn(id,3,hit.pos.y);
            ana->FillNtupleFColumn(id,4,hit.pos.z);
            ana->FillNtupleFColumn(id,5,hit.time);
            ana->FillNtupleFColumn(id,6,hit.mom.x);
            ana->FillNtupleFColumn(id,7,hit.mom.y);
            ana->FillNtupleFColumn(id,8,hit.mom.z);
            ana->FillNtupleFColumn(id,9,hit.pol.x);
            ana->FillNtupleFColumn(id,10,hit.pol.y);
            ana->FillNtupleFColumn(id,11,hit.pol.z);
            ana->FillNtupleFColumn(id,12,hit.wavelength);
            ana->AddNtupleRow(id);
        }
#endif
    }

  void SensorSD::GEANT4Hits(G4Step* step){
      auto analysisManager = G4AnalysisManager::Instance();
      G4int id;
      auto track=step->GetTrack();
      auto event=G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
      if(this->GetName()=="Camera")  id=0;
      else id=4;
      analysisManager->FillNtupleDColumn(id,0, event);
      analysisManager->FillNtupleDColumn(id,1, track->GetTrackID());
      analysisManager->FillNtupleDColumn(id,2, track->GetProperTime()/CLHEP::ns);
      analysisManager->FillNtupleDColumn(id,3, track->GetPosition()[0]/CLHEP::mm);
      analysisManager->FillNtupleDColumn(id,4, track->GetPosition()[1]/CLHEP::mm);
      analysisManager->FillNtupleDColumn(id,5, track->GetPosition()[2]/CLHEP::mm);
      analysisManager->AddNtupleRow(id);

  }
  void SensorSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
  {
    //  int HCID = G4SDManager::GetSDMpointer()->
    //    GetCollectionID(this->GetCollectionName(0));
    //  // }
    // HCE->AddHitsCollection(HCID, HC_);


  }


} // end namespace nexus
