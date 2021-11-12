#ifndef EventAction_hh
#define EventAction_hh 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"
#include <vector>


class G4VPhysicalVolume;
class SteppingAction;
class G4Event;


class EventAction : public G4UserEventAction {
 public:
  EventAction();
  EventAction(DetectorConstruction* dc);
  ~EventAction();

 public:
  void BeginOfEventAction(const G4Event *);
  void EndOfEventAction(const G4Event *);
  

 private:
  DetectorConstruction* fdetCon;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
