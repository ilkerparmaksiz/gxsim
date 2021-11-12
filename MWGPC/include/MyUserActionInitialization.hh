#ifndef MyUserActionInitialization_hh
#define MyUserActionInitialization_hh

#include "G4VUserActionInitialization.hh"
#include "DetectorConstruction.hh"

//class DetectorConstruction;
class PhysicsList;

class MyUserActionInitialization : public G4VUserActionInitialization{
	public:
	MyUserActionInitialization();
        MyUserActionInitialization(DetectorConstruction *);
        ~MyUserActionInitialization();
	
	void Build() const;
	void BuildForMaster() const;
	
	private:
        DetectorConstruction *fdetCon;
};

#endif
