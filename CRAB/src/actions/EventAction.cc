#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "RunAction.hh"
#include "config.h"

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "Analysis.hh"
#include "SteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "DegradModel.hh"
#include "GarfieldVUVPhotonModel.hh"
#ifdef With_Opticks
#include "SensorSD.hh"
#  include "SEvt.hh"
#  include "NP.hh"
#  include "G4CXOpticks.hh"
namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
using namespace sensorsd;
#endif
EventAction::EventAction() {
  
}

EventAction::~EventAction() {
	G4cout << "Deleting EventAction" << G4endl;
}


void EventAction::BeginOfEventAction(const G4Event *ev) {
  G4cout << " EventAction::BeginOfEventAction() "<< ev->GetEventID() << G4endl;
    DegradModel* dm = (DegradModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("DegradModel"));
    if(dm)
        dm->Reset();

    G4PrimaryVertex* pVtx;
    pVtx = ev->GetPrimaryVertex();
    if (pVtx)
      {
        G4double PKE = pVtx->GetPrimary(0)->GetKineticEnergy();
	dm->SetPrimaryKE(PKE);
      }

    fEDepPrim = 0.0;
}

void EventAction::EndOfEventAction(const G4Event *evt) {

    GarfieldVUVPhotonModel* gvm = (GarfieldVUVPhotonModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("GarfieldVUVPhotonModel"));
    if(gvm)
      gvm->Reset(); // zero out the sensor: meaning reset the nexcitations, which is cumulative.
    G4cout << " EventAction::EndOfEventAction()  "  << evt->GetEventID()<< G4endl;
#ifdef With_Opticks

    G4cout<<" Opticks End of Event Action" <<G4endl;
    G4AutoLock lock(&opticks_mt);
    G4CXOpticks * g4cx=G4CXOpticks::Get();

    G4int eventID=evt->GetEventID();
    G4int ngenstep=SEvt::GetNumGenstepFromGenstep(eventID);
    G4int nphotons=SEvt::GetNumPhotonCollected(eventID);



    // Simulate the photons
      if(nphotons>0 and ngenstep>0){
          std::cout<<g4cx->desc()<<std::endl;
          std::cout<<"--- G4Optickx ---" << g4cx->descSimulate() <<std::endl;
          g4cx->simulate(eventID,0); // For Simulation
          //cudaDeviceSynchronize();
          //g4cx->render();  // For Rendering

      }


    //SensorSD* PMT = (SensorSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector("PMT");
    SensorSD* Camera = (SensorSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector("/Sensor/Camera");

    G4cout << "Number of Steps Generated " <<ngenstep << G4endl;
    G4cout << "Number of Photons Generated " <<nphotons << G4endl;
    G4cout << "Number of Hits Opticks  " <<SEvt::GetNumHit(eventID)<< G4endl;

    if(SEvt::GetNumHit(eventID)>0){
        //PMT->OpticksHits();
        Camera->OpticksHits();
        G4CXOpticks::Get()->reset(eventID);

    }


#endif

    G4int id(3);
    G4double PPID = 0.; G4double PKE = 0.;
    G4PrimaryVertex* pVtx;
    pVtx = evt->GetPrimaryVertex();
    if (pVtx)
      {
	PKE = pVtx->GetPrimary(0)->GetKineticEnergy();
	PPID = pVtx->GetPrimary(0)->GetPDGcode();
      }

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int  event = evt->GetEventID();
    G4int row(0);
    analysisManager->FillNtupleDColumn(id,row, event); row++;
    analysisManager->FillNtupleDColumn(id,row, (G4double)PPID); row++;
    analysisManager->FillNtupleDColumn(id,row, PKE); row++;
    analysisManager->FillNtupleDColumn(id,row, fEDepPrim); row++;
    analysisManager->AddNtupleRow(id);



}

void EventAction::EDepPrim(const G4double &Ed)
{
  fEDepPrim+=Ed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
