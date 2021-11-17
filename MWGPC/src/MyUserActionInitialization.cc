#include "MyUserActionInitialization.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "GasBoxSD.hh"
#include "SteppingAction.hh"

MyUserActionInitialization::MyUserActionInitialization(){}
MyUserActionInitialization::MyUserActionInitialization(DetectorConstruction *det){ fdetCon = det; std::cout << "MyUserActionInitialization::Constructor():  " << fdetCon<< std::endl;}

MyUserActionInitialization::~MyUserActionInitialization(){}

void MyUserActionInitialization::Build() const {
	PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
	SetUserAction(primary);
	SteppingAction* stepAct = new SteppingAction();
	SetUserAction(stepAct);

	std::cout << "MyUserActionInitialization::Build(): Registering EventAction() " << fdetCon << std::endl;
	EventAction* evt = new EventAction(fdetCon);
	SetUserAction(evt);
	std::cout << "MyUserActionInitialization::Build(): Registering RunAction() " << fdetCon<< std::endl;	
	SetUserAction(new RunAction(fdetCon));
}

void MyUserActionInitialization::BuildForMaster() const {
	SetUserAction(new RunAction(fdetCon));
}
