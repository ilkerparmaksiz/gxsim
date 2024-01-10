#include "GasBoxSD.hh"
#include "G4Region.hh"
#include "G4String.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "DetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalPhoton.hh"
#include "S2Photon.hh"
#include "NEST/G4/NESTProc.hh"
GasBoxSD::GasBoxSD(G4String name) : G4VSensitiveDetector(name),
    fXenonHitsCollection(NULL), fGarfieldExcitationHitsCollection(NULL){
    collectionName.insert("XHC");
    collectionName.insert("GEHC");
    
    XHCID=-1;
    GEHCID=-1;
}

GasBoxSD::~GasBoxSD(){}


void GasBoxSD::Initialize(G4HCofThisEvent * HCE){
    fXenonHitsCollection = new XenonHitsCollection(SensitiveDetectorName, collectionName[0]);
    fGarfieldExcitationHitsCollection = new GarfieldExcitationHitsCollection(SensitiveDetectorName, collectionName[1]);
    if(XHCID==-1){
        XHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        GEHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]);
    }
    HCE->AddHitsCollection(XHCID,fXenonHitsCollection);
    HCE->AddHitsCollection(GEHCID,fGarfieldExcitationHitsCollection);

    G4cout << "!!!!!!!!!!! GasBoxSD Intialized  !!!!!!!!!!!!!" << G4endl;
}

G4bool GasBoxSD::ProcessHits(G4Step* aStep, G4TouchableHistory* hist){
    G4Track* aTrack = aStep->GetTrack();
    G4ParticleDefinition *currDeff =aTrack->GetDefinition();
    G4double edep=aStep->GetTotalEnergyDeposit();
    // do not worry if there is no energy in the detector
    if (edep<=0) return false;
    //
    if (currDeff==G4OpticalPhoton::Definition() || currDeff==S2Photon::Definition() || currDeff==NEST::NESTThermalElectron::Definition()) return false;
    G4bool isPrint=false;

    if(isPrint){
        G4cout<<"Processing the hit after" <<G4endl;
        G4cout << "GasBox Hit!!" << G4endl;
        G4cout << "Particle ID: " << aTrack->GetTrackID() << G4endl;
        G4cout << "Particle Name: " << currDeff->GetParticleName() << G4endl;
        G4cout << "Particle electron: " << aTrack->GetKineticEnergy() << G4endl;

        G4cout << "Particle Deposit: " << edep << G4endl;
    }
    return true;

    
}

void GasBoxSD::EndOfEvent (G4HCofThisEvent * hce){
  
  auto HC = static_cast<XenonHitsCollection*>(hce->GetHC(XHCID));
    int entries = HC->entries();
    G4cout << "Number of Electrons: " << entries << G4endl;
  /*
    for(int i=0;i<entries;i++){
        auto hit = (*HC)[i];
        G4cout << hit->GetPos() << " " << hit->GetTime() << G4endl;
    }
    */
  /*
    auto HC1 = static_cast<GarfieldExcitationHitsCollection*>(hce->GetHC(GEHCID));
    int entries1 = HC1->entries();
    G4cout << "GasBoxSD::EndOfEvent (): Number of De-excitations: " << entries1 << G4endl;
  */
    /*
    for(int i=0;i<entries1;i++){
        auto hit = (*HC1)[i];
        G4cout << hit->GetPos() << " " << hit->GetTime() << G4endl;
    }
    */
    DrawAll();
}

void GasBoxSD::DrawAll(){}
