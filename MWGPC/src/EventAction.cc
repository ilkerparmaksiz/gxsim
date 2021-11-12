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

  // get the DetectorConstruction Heed DeltaElectron Tracker class instance.
  auto hDEtracks = fdetCon->fHDEt; // prob should create/use a getter .... EC, 11-Nov2022.
  
  for (int bin = 0; bin<hDEtracks->GetMaxbin(); bin++)
      {
        if (hDEtracks->GetSensor()->GetSignal("s2", bin) != 0.)
          std::cout << " wire electron signal: " << bin*hDEtracks->GetBinsz() << " nsec: "<< hDEtracks->GetSensor()->GetSignal("s2", bin) << std::endl;
	analysisManager->FillP1(0,(bin+bin+1)/2.*hDEtracks->GetBinsz(),hDEtracks->GetSensor()->GetSignal("s2", bin));
      }


  hDEtracks->GetSensor()->ClearSignal();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
