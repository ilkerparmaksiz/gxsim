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

  auto analysisManager = G4AnalysisManager::Instance();

  G4int evID (ev->GetEventID());
  std::string p1name("WireSig_");
  analysisManager->CreateP1(p1name+std::to_string(evID),"ULBPC Wire Signal [fC/nsec]", fdetCon->fNumBins, 0., fdetCon->fNumBins*fdetCon->fBinSz,-1E5,+1E5);  // profile histo in bins of 2000/100 nsec

}

void EventAction::EndOfEventAction(const G4Event *evt) {
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
