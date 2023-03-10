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
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4ProcessTable.hh"
#include "G4RadioactiveDecay.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
  : G4VPrimaryGenerator(), momentum_{}
{

  msg_ = new G4GenericMessenger(this, "/Generator/SingleParticle/",
    "Control commands of single-particle generator.");

  msg_->DeclarePropertyWithUnit("momentum", "mm",  momentum_, "Set particle 3-momentum.");

  //  --- Get Xenon file --- 
	char* nexus_path = std::getenv("CRABPATH");
	if (nexus_path == nullptr) {
		G4Exception("[PrimaryGenerator]", "Constructor()", FatalException,
					"Environment variable CRABPATH not defined!");
	}

	std::string crab_path(nexus_path);



  FileHandler.GetEvent(crab_path + "/data/pb210_electron.csv", electron_data);
  FileHandler.GetEvent(crab_path + "/data/pb210_alpha.csv",    alpha_data);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertexOpt(G4Event* event, std::vector<double> &xyzb)
{
  //vertex A uniform on a cylinder
  //
  G4int n_particle = 1;

  G4double KE = 1 * MeV;

  G4ThreeVector positionA( xyzb.at(0), xyzb.at(1), xyzb.at(2));
  G4double timeA = 0*s;

  G4ThreeVector e_momentum = {};
  G4ThreeVector alpha_momentum = {};

  G4double e_Ke = 1.16*MeV;
  G4double alpha_Ke = 5*MeV;

  G4bool useNeedle = true; 

  // Generate events off the surface of the needle
  if (useNeedle){

    
    G4double maxRad_ = (0.42)*mm + 2*nm;
    G4double halfLength = 1 * mm;
    G4double iniPhi_ = 0;
    G4double deltaPhi_ = twopi;
    G4ThreeVector origin_ = {-1.6*cm - 1*mm, 0, - 5*cm };
    
    G4RotationMatrix* rotateHolder = new G4RotationMatrix();
    rotateHolder->rotateY(90.*deg);


    G4double phi = (iniPhi_ + (G4UniformRand() * deltaPhi_));
    G4double rad = maxRad_;

    positionA = {rad * cos(phi), rad * sin(phi), (G4UniformRand() * 2.0 - 1.0) * halfLength  };
  
    positionA *= *rotateHolder;
    
    // Translating
    positionA += origin_;
    
  }

  G4bool useEventsFromFile = true; 

  if (useEventsFromFile){

    // Get some random index
    G4int event  = round(G4UniformRand()*electron_data.size());
    std::vector<G4double> electron_event = electron_data[event];
    std::vector<G4double> alpha_event = alpha_data[event];

    e_momentum = {electron_event[0], electron_event[1], electron_event[2]};
    e_momentum = e_momentum.unit();

    alpha_momentum = {alpha_event[0], alpha_event[1], alpha_event[2]};
    alpha_momentum = alpha_momentum.unit();

    // std::cout << "Sampled electron with px, py, pz, KE:" << electron_event[0] << ", " << electron_event[1] << ", " << electron_event[2] << ", " << electron_event[3] << std::endl;
    // std::cout << "Sampled alpha with px, py, pz, KE:" << alpha_event[0] << ", " << alpha_event[1] << ", " << alpha_event[2] << ", " << alpha_event[3] << std::endl;

    std::cout << "Sampled electron with px, py, pz, KE:" << e_momentum.x() << ", " << e_momentum.y()  << ", " << e_momentum.z()  << ", " << electron_event[3] << std::endl;
    std::cout << "Sampled alpha with px, py, pz, KE:" <<alpha_momentum.x()<< ", " << alpha_momentum.y() << ", " << alpha_momentum.z() << ", " << alpha_event[3] << std::endl;

    e_Ke     = electron_event[3]*MeV;
    alpha_Ke = alpha_event[3]*MeV;

  }
  else {
    e_momentum     = momentum_.unit();
    alpha_momentum = momentum_.unit();

  }


  // 
  G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);

  G4ParticleDefinition* particleDefinition;
  G4PrimaryParticle* particle1;
  G4double mass;
  G4double energy;
  G4double pmod;
  G4ThreeVector p;


  for (int ii = 0; ii < n_particle; ii++) {
      
    // Particle 1 at vertex A
    
    // Initialise the alpha
    if (ii == 0){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
      particle1 = new G4PrimaryParticle(particleDefinition);
      KE = alpha_Ke;

      mass   = particleDefinition->GetPDGMass();
      energy = KE + mass;
      pmod = std::sqrt(energy*energy - mass*mass);
      p = pmod * alpha_momentum;
      particle1->SetMomentum(p.x(), p.y(), p.z());
      std::cout << "\nPrimaryGenerator: Adding alpha with " << particle1->GetKineticEnergy()/keV << " keV w ux,uy,uz " << p.x() << ", " << p.y() << ", " << p.z()<< " to vertexA."  << std::endl;
     
    }
    // Initalise the electron
    else {
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
      particle1 = new G4PrimaryParticle(particleDefinition);
      KE = e_Ke;

      mass   = particleDefinition->GetPDGMass();
      energy = KE + mass;
      pmod = std::sqrt(energy*energy - mass*mass);
      p = pmod * e_momentum;
      particle1->SetMomentum(p.x(), p.y(), p.z());
      std::cout << "PrimaryGenerator: Adding electron with " << particle1->GetKineticEnergy()/keV << " keV  e- w ux,uy,uz " << p.x() << ", " << p.y() << ", " << p.z()<< " to vertexA.\n"  << std::endl;
    }
    
    // Add particle to the vertex
    vertexA->SetPrimary(particle1);
    
  }

  event->AddPrimaryVertex(vertexA);
  //  SetParticlePosition(positionA);
  //  std::cout << "PrimaryGenerator: Added " << n_particle << " isotropic electrons as primaries at " << x/1000. <<", " << y/1000. << ", " << z/1000. << " [m]. "  << std::endl;
  vertexA->Print();
}

// -------------------


void PrimaryGenerator::GeneratePrimaryVertexIon(G4Event* event, std::vector<double> &xyzb){


  G4ThreeVector positionA( xyzb.at(0), xyzb.at(1), xyzb.at(2));
  G4double timeA = 0*s;

  G4bool useNeedle = true; 

  // Generate events off the surface of the needle
  if (useNeedle){

    
    G4double maxRad_ = (0.42)*mm + 2*nm;
    G4double halfLength = 1 * mm;
    G4double iniPhi_ = 0;
    G4double deltaPhi_ = twopi;
    G4ThreeVector origin_ = {-1.6*cm - 1*mm, 0, - 5*cm };
    
    G4RotationMatrix* rotateHolder = new G4RotationMatrix();
    rotateHolder->rotateY(90.*deg);


    G4double phi = (iniPhi_ + (G4UniformRand() * deltaPhi_));
    G4double rad = maxRad_;

    positionA = {rad * cos(phi), rad * sin(phi), (G4UniformRand() * 2.0 - 1.0) * halfLength  };
  
    positionA *= *rotateHolder;
    
    // Translating
    positionA += origin_;
    
  }

  G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);

  G4ParticleDefinition* particleDefinition;
  G4PrimaryParticle* particle1;
  G4double mass;
  G4double energy;
  G4double pmod;
  G4ThreeVector p;

  // Initialise the alpha
  particleDefinition = G4IonTable::GetIonTable()->GetIon(82, 210, 0.);
  particleDefinition->SetPDGLifeTime(1.*ps); 
  G4PrimaryParticle* ion = new G4PrimaryParticle(particleDefinition);
  
  // Add particle to the vertex
  vertexA->SetPrimary(ion);
  event->AddPrimaryVertex(vertexA);
  vertexA->Print();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


