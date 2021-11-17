#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "RunAction.hh"

#include "G4SDManager.hh"
#include "G4Threading.hh"

#include "Analysis.hh"
#include "SteppingAction.hh"

EventAction::EventAction() {
  
}

EventAction::EventAction(DetectorConstruction* dc) {fdetCon = dc;};

EventAction::~EventAction() {
	G4cout << "Deleting EventAction" << G4endl;
}


void EventAction::BeginOfEventAction(const G4Event *ev) {

}

void EventAction::EndOfEventAction(const G4Event *evt) {

  auto analysisManager = G4AnalysisManager::Instance();
  //  analysisManager->ChangeDirectory("histo");

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
