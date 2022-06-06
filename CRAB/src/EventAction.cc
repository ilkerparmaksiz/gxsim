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
    G4cout << " EventAction::BeginOfEventAction()  1 " << G4endl;
}

void EventAction::EndOfEventAction(const G4Event *evt) {

     G4cout << " EventAction::EndOfEventAction()  " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
