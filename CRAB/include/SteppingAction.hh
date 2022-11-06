#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4OpBoundaryProcess.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction(EventAction *eva) : fEventAction(eva) {};
  ~SteppingAction(){};

  void UserSteppingAction(const G4Step *);
 
 private:
  EventAction* fEventAction;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
