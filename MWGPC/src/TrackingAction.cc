//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
//#include "Run.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4HadronicProcessType.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* evt)
 :G4UserTrackingAction(), fDetector(det), fEventAction(evt)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* g4t)
{
  
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int evID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4int pntid = g4t->GetParentID();
  G4int tid = g4t->GetTrackID();

  if (tid == 1 && pntid == 0)
    {
      G4int id(0);
      analysisManager->FillNtupleDColumn(id,0, evID);
      analysisManager->FillNtupleDColumn(id,1, g4t->GetParticleDefinition()->GetPDGEncoding());
      analysisManager->FillNtupleDColumn(id,2, g4t->GetKineticEnergy()/keV);
      analysisManager->FillNtupleDColumn(id,3, g4t->GetPosition()[0]/mm);
      analysisManager->FillNtupleDColumn(id,4, g4t->GetPosition()[1]/mm);
      analysisManager->FillNtupleDColumn(id,5, g4t->GetPosition()[2]/mm);
      analysisManager->FillNtupleDColumn(id,6, g4t->GetMomentumDirection()[0]);
      analysisManager->FillNtupleDColumn(id,7, g4t->GetMomentumDirection()[1]);
      analysisManager->FillNtupleDColumn(id,8, g4t->GetMomentumDirection()[2]);
    }

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* )
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

