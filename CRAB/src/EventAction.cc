#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "RunAction.hh"

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "SteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "DegradModel.hh"
#include "GarfieldVUVPhotonModel.hh"

EventAction::EventAction() {
  
}

EventAction::~EventAction() {
	G4cout << "Deleting EventAction" << G4endl;
}


void EventAction::BeginOfEventAction(const G4Event *ev) {
  G4cout << " EventAction::BeginOfEventAction()  0 " << G4endl;
    DegradModel* dm = (DegradModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("DegradModel"));
    if(dm)
        dm->Reset();

    fEDepPrim = 0.0;
    G4cout << " EventAction::BeginOfEventAction()  1 " << G4endl;
}

void EventAction::EndOfEventAction(const G4Event *evt) {

    GarfieldVUVPhotonModel* gvm = (GarfieldVUVPhotonModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("GarfieldVUVPhotonModel"));
    if(gvm)
      gvm->Reset(); // zero out the sensor: meaning reset the nexcitations, which is cumulative.

    G4cout << " EventAction::EndOfEventAction()  " << G4endl;


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
