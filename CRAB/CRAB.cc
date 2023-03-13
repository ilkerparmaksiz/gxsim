/**
 *\file CRAB.cc
 *\brief Main program of CRAB simulation
 *\author Claire Couratin
 *\date December 2014


You can find the Wiki on the subversion repository:
https://svs.icts.kuleuven.be/projects/svs_project014/wiki/Wiki
 */
#include <ctime>
#include <cstdio>
#include <time.h>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VSteppingVerbose.hh"
#include "Randomize.hh" 

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "MyUserActionInitialization.hh"
#include "GasModelParameters.hh"

int main(int argc, char** argv) {
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager();
  //    runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
  runManager->SetNumberOfThreads(1);
#else
    G4RunManager* runManager = new G4RunManager();
    G4cout << "Creation of G4RunManager" << G4endl;
#endif


  
  
  G4int randseed = atoi(argv[2]);
  G4Random::setTheSeed(randseed);
  G4cout << "Setting the Random seed: " << randseed << G4endl;
  
  G4cout << "Creation of the gas model parameter class" << G4endl;
  GasModelParameters* gmp = new GasModelParameters();
    
  G4cout << "Creation of DetectorConstruction" << G4endl;
  DetectorConstruction* detector = new DetectorConstruction(gmp);
  runManager->SetUserInitialization(detector);

  
  G4cout << "Creation of PhysicsList" << G4endl;
  PhysicsList* physics = new PhysicsList();
  runManager->SetUserInitialization(physics);
  
  runManager->SetUserInitialization(new MyUserActionInitialization());
 
  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive* ui = nullptr;
  G4VisManager* visManager = nullptr;
  //runManager->Initialize();

  if (argc == 1)  //! define UI terminal for interactive mode:
  {
    ui = new G4UIExecutive(argc, argv);
    visManager = new G4VisExecutive();
    visManager->Initialize();
    //#ifdef G4UI_USE

    //#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
    //#endif

    ui->SessionStart();
    delete ui;
    //#endif
    //#ifdef G4VIS_USE

    //#endif

  } else  //! batch mode:
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
  
    if (argc < 3) {
      G4cout << "No random seed has been provided" << G4endl;
      delete runManager;
      return 0;
    }
    G4cout << "About to launch: " << command << " " << fileName << G4endl;
    UImanager->ApplyCommand(command + fileName);

    time_t start=time(0);

    double duration = difftime(time(0),start);
    cout << "Simulation Time: " << duration << endl;
  }

  if (visManager)
    delete visManager;
  
  // if (runManager)
  //   {
  //     std::cout << "Deleting runManager" << std::endl;
  //     delete runManager;
  //   }
  
   return 0;
}
