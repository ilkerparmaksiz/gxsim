#include "SteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "DetectorConstruction.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "GasBoxSD.hh"



SteppingAction::SteppingAction(){

}


void SteppingAction::UserSteppingAction(const G4Step *aStep) {


G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
 G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
G4TouchableHistory* theTouchable
 = (G4TouchableHistory*)(preStepPoint->GetTouchable());
G4ThreeVector worldPrePos = preStepPoint->GetPosition();
 G4ThreeVector worldPostPos = postStepPoint->GetPosition();

 // std::cout << "UserSteppingAction() preStep worldPos,volName is "  << worldPrePos << ", " << preStepPoint-> GetPhysicalVolume () -> GetName ()  << std::endl;
 // std::cout << "UserSteppingAction() postStep worldPos,volName is " <<  worldPostPos<< ", " << postStepPoint-> GetPhysicalVolume () -> GetName ()  << std::endl;
 
}
