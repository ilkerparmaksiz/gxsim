//
//  PrimaryGeneratorAction.cpp
//  Xenon
//
//  Created by Lennert De Keukeleere on 25/10/2018.
//

#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4ThreeVector.hh"
#include "G4Geantino.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(){
    fparticleGun = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fparticleGun;
    G4cout << "Deleting PrimaryGeneratorAction" << G4endl;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  std::vector<double> xyzbounds { 
      fparticleGun->GetCurrentSource()->GetPosDist()->GetHalfX(),
      fparticleGun->GetCurrentSource()->GetPosDist()->GetHalfY(),
      fparticleGun->GetCurrentSource()->GetPosDist()->GetHalfZ()};

  std::cout << "GeneratePrimaries() particle requested is " << fparticleGun->GetParticleDefinition() << std::endl;
  std::cout << "GeneratePrimaries() vtx bounds:  " << xyzbounds[0] << ", " << xyzbounds[1] << ", " << xyzbounds[2] << std::endl;


  fparticleGun->GeneratePrimaryVertex(anEvent);

}

