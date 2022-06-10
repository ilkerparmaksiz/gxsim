//Runaction
#include "RunAction.hh"
#include "G4Run.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "EventAction.hh"

#include "GasBoxSD.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

RunAction::RunAction(){
  G4cout << "Creating AnalysisManager" << G4endl;
  auto analysisManager = G4AnalysisManager::Instance();
//  analysisManager->SetNtupleMerging(true,0,0,10000000);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  
  analysisManager->SetFileName("output.root"); 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");
  
  analysisManager->SetNtupleActivation(true);

  G4cout << "Creating RunAction" << G4endl;
}

RunAction::~RunAction() { 
	G4cout << "Deleting RunAction" << G4endl;
	delete G4AnalysisManager::Instance();  
}


void RunAction::BeginOfRunAction(const G4Run* aRun) {
  G4Random::showEngineStatus();

  G4cout << "Starting run " << aRun->GetRunID() << G4endl;
  time_t currentTime;
  tm* ptm;
  time(&currentTime);
  ptm = localtime(&currentTime);
  G4cout << "Time: " << asctime(ptm) << G4endl;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
  
  analysisManager->CreateNtuple("PMT", "PMT Hits"); 
  analysisManager->CreateNtupleDColumn("Event");         //column 1
  analysisManager->CreateNtupleDColumn("PID");         //column 2
  analysisManager->CreateNtupleDColumn("Time");         //column 3
  analysisManager->CreateNtupleDColumn("X");         //column 4
  analysisManager->CreateNtupleDColumn("Y");         //column 5
  analysisManager->CreateNtupleDColumn("Z");         //column 6
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("S1", "S1 particles"); 
  analysisManager->CreateNtupleDColumn("Event");         //column 1
  analysisManager->CreateNtupleDColumn("PID");         //column 2
  analysisManager->CreateNtupleDColumn("Time");         //column 3
  analysisManager->CreateNtupleDColumn("X");         //column 4
  analysisManager->CreateNtupleDColumn("Y");         //column 5
  analysisManager->CreateNtupleDColumn("Z");         //column 6
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("S2", "S2 particles"); 
  analysisManager->CreateNtupleDColumn("Event");         //column 1
  analysisManager->CreateNtupleDColumn("PID");         //column 2
  analysisManager->CreateNtupleDColumn("Time");         //column 3
  analysisManager->CreateNtupleDColumn("X");         //column 4
  analysisManager->CreateNtupleDColumn("Y");         //column 5
  analysisManager->CreateNtupleDColumn("Z");         //column 6
  analysisManager->FinishNtuple();

  analysisManager->SetNtupleActivation(true);
}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "End of run OK!" << G4endl;
  time_t currentTime;
  tm* ptm;
  time(&currentTime);
  ptm = localtime(&currentTime);
  G4cout << "Simulation finished." << G4endl << "Time: " << asctime(ptm)
         << G4endl;
  G4Random::showEngineStatus();
}
