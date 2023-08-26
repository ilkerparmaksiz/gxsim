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
#include "FileHandling.hh"
#include "G4OpticalPhoton.hh"
#include "S2Photon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
  : G4VPrimaryGenerator(),momentum_(1,1,1), ParticleType_("opticalphoton"),Position_(0),energy_(0),Iso_(true),useNeedle(true),fAmount_(1),GeneratorMode_("Single")
{

  msg_ = new G4GenericMessenger(this, "/Generator/SingleParticle/",
    "Control commands of single-particle generator.");

  msg_->DeclareProperty("momentum",  momentum_, "Set particle 3-momentum.");
  msg_->DeclareProperty("energy",  energy_, "Set particle energy");
  msg_->DeclareProperty("ParticleType",  ParticleType_, "Set particle type alpha,e- or etc");
  msg_->DeclarePropertyWithUnit("pos", "cm",  Position_, "Set Position x,y,z");
  msg_->DeclareProperty("Isotropic",  Iso_, "Isotropic Distribution");
  msg_->DeclareProperty("useNeedle",  useNeedle, "Isotropic Distribution");
  msg_->DeclareProperty("Mode",  GeneratorMode_, "The mode of the generator to run");
  msg_->DeclareProperty("Amount",  fAmount_, "Amount of particles at a time");

    //msg_->DeclareProperty("pos", "cm",  Position_, "Set Position x,y,z");
  //  --- Get Xenon file --- 
	char* nexus_path = std::getenv("CRABPATH");
	if (nexus_path == nullptr) {
		G4Exception("[PrimaryGenerator]", "Constructor()", FatalException,
					"Environment variable CRABPATH not defined!");
	}

    CRABPATH= getenv("CRABPATH");

  //FileHandler.GetEvent(crab_path + "/data/pb210_electron.csv", electron_data);
  //FileHandler.GetEvent(crab_path + "/data/pb210_alpha.csv",    alpha_data);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

// -------------------
void PrimaryGenerator::Generate(G4Event* event, std::vector<double> &xyz){
  // Set Needle Path
  using namespace filehandler;
  FileHandling *f=new FileHandling();
  NeedlePointPath= CRABPATH +"data/"+GeneratorMode_+".txt";

  if (GeneratorMode_ == "Single"){
    std::cout <<"Generating with Single Particle Mode for event: " << event->GetEventID() << std::endl;
      if (fAmount_==0) fAmount_=1;
      for (int i=0; i<fAmount_;i++) GenerateSingleParticle(event);
  }
  else if(GeneratorMode_!=" " && f->FileCheck(NeedlePointPath)){
      NeedlePoints=f->GetThreeVectorData(NeedlePointPath,',',1);
      if (fAmount_==0) fAmount_=1;
      for (int i=0; i<fAmount_;i++) GenerateFromSurface(event);
  }
  else {
    std::cout <<"Generating with Ion Mode for event: " << event->GetEventID() << std::endl;
    GeneratePrimaryVertexIon(event, xyz);
  }

}



// -------------------

void PrimaryGenerator::GeneratePrimaryVertexIon(G4Event* event, std::vector<double> &xyzb){


  G4ThreeVector positionA( xyzb.at(0), xyzb.at(1), xyzb.at(2));
  G4double timeA = 0*s;

  // Generate events off the surface of the needle
  if (useNeedle){

    
    G4double maxRad_ = (0.42)*mm + 2*nm;
    G4double halfLength = 1 * mm;
    G4double iniPhi_ = 0;
    G4double deltaPhi_ = twopi;
    G4ThreeVector origin_ = {0 , -1.6*cm-1*mm, - 5.25*cm };
    
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

  pmod = 500*keV;
  p = pmod * momentum_;

  // Initialise the alpha
  particleDefinition = G4IonTable::GetIonTable()->GetIon(84, 210, 0.); // Po210 decay of 5.3 MeV alpha
  particleDefinition->SetPDGLifeTime(1.*ps); 
  G4PrimaryParticle* ion = new G4PrimaryParticle(particleDefinition);
  ion->SetMomentum(p.x(), p.y(), p.z());
  
  // Add particle to the vertex
  vertexA->SetPrimary(ion);
  event->AddPrimaryVertex(vertexA);
  vertexA->Print();
}

void PrimaryGenerator::GenerateSingleParticle(G4Event * event) {


    G4double NeedleOffset=1*mm;
    G4ThreeVector positionA( Position_[0], Position_[1], Position_[2]);
    energy_=energy_*MeV;
    G4ParticleDefinition* particleDefinition;
    G4PrimaryParticle* particle1;
    G4double mass;
    G4double energy;
    G4double pmod;
    G4ThreeVector p;



    // Particle 1 at vertex A

    // Initialise the Single Particle

    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ParticleType_);
    particle1 = new G4PrimaryParticle(particleDefinition);

    mass   = particleDefinition->GetPDGMass();
    energy = energy_ + mass;
    pmod = std::sqrt(energy*energy - mass*mass);

    if(Iso_){
        G4ThreeVector SphericalCoord;
        G4double iniPhi=0;
        G4double deltaPhi_=pi;
        G4double deltatheta=twopi;

        G4double phi = (iniPhi + (G4UniformRand() * deltaPhi_));
        G4double theta=(iniPhi+(G4UniformRand() * deltatheta));
        G4double rad = 1;
        SphericalCoord ={rad*sin(phi)*cos(theta),rad*sin(phi)*sin(theta),rad*cos(phi)};
        p=SphericalCoord*pmod;

    }else{
        p = pmod * momentum_;

    }

    particle1->SetMomentum(p.x(), p.y(), p.z());
    //std::cout << "\nPrimaryGenerator: Adding particle with " << particle1->GetKineticEnergy()/keV << " keV w ux,uy,uz " << p.x() << ", " << p.y() << ", " << p.z()<< " to vertexA."  << std::endl;

    // Generate events off the surface of the needle
    if (useNeedle){

        // Shift y position
        positionA = {Position_[0], Position_[1]- NeedleOffset, Position_[2]};


        G4double maxRad_ = (0.42)*mm + 2*nm;
        G4double halfLength = 1 * mm;
        G4double iniPhi_ = 0;
        G4double deltaPhi_ = twopi;

        G4ThreeVector origin_ = {0 , -1.6*cm-1*mm, - 5.25*cm };

        G4RotationMatrix* rotateHolder = new G4RotationMatrix();
        rotateHolder->rotateY(90.*deg);


        G4double phi = (iniPhi_ + (G4UniformRand() * deltaPhi_));
        G4double rad = maxRad_;

        positionA = {rad * cos(phi), rad * sin(phi), (G4UniformRand() * 2.0 - 1.0) * halfLength  };

        positionA *= *rotateHolder;

        // Translating
        positionA += origin_;

    }

    G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA,0);


    RandomPolarization(particle1);

    // Add particle to the vertex
    vertexA->SetPrimary(particle1);
    event->AddPrimaryVertex(vertexA);
    //vertexA->Print();
}
void PrimaryGenerator::GenerateFromSurface(G4Event* evt){
    G4int Size=NeedlePoints.size()-1;
    G4int Index;
    Index= round(G4UniformRand()*Size);

    G4double NeedleOffset=1*mm;
    G4ThreeVector positionA=NeedlePoints[Index];
    energy_=energy_*MeV;
    G4ParticleDefinition* particleDefinition;
    G4PrimaryParticle* particle1;
    G4double mass;
    G4double energy;
    G4double pmod;
    G4ThreeVector p;



    // Particle 1 at vertex A

    // Initialise the Single Particle

    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ParticleType_);
    particle1 = new G4PrimaryParticle(particleDefinition);

    mass   = particleDefinition->GetPDGMass();

    energy = energy_ + mass;
    pmod = std::sqrt(energy*energy - mass*mass);

    if(Iso_){
        G4ThreeVector SphericalCoord;
        G4double iniPhi=0;
        G4double deltaPhi_=pi;
        G4double deltatheta=twopi;

        G4double phi = (iniPhi + (G4UniformRand() * deltaPhi_));
        G4double theta=(iniPhi+(G4UniformRand() * deltatheta));
        G4double rad = 1;
        SphericalCoord ={rad*sin(phi)*cos(theta),rad*sin(phi)*sin(theta),rad*cos(phi)};
        p=SphericalCoord*pmod;

    }else{
        p = pmod * momentum_;

    }

    particle1->SetMomentum(p.x(), p.y(), p.z());
    //std::cout << "\nPrimaryGenerator: Adding particle with " << particle1->GetKineticEnergy()/keV << " keV w ux,uy,uz " << p.x() << ", " << p.y() << ", " << p.z()<< " to vertexA."  << std::endl;
    //G4cout<<"Particle Position is " << positionA <<G4endl;
    G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA,0);
    // This is needed for optical photons
    RandomPolarization(particle1);
    // Add particle to the vertex
    vertexA->SetPrimary(particle1);
    evt->AddPrimaryVertex(vertexA);


}
// This is to create polarization for optical photons and S2Photons
void PrimaryGenerator::RandomPolarization(G4PrimaryParticle *particle) {
    if(particle->GetParticleDefinition()==G4OpticalPhoton::Definition() or particle->GetParticleDefinition()==S2Photon::Definition() ) {
        G4double angle = G4UniformRand() * 360.0 * deg;
        G4ThreeVector normal(1., 0., 0.);
        G4ThreeVector kphoton = particle->GetMomentumDirection();
        G4ThreeVector product = normal.cross(kphoton);
        G4double modul2 = product * product;

        G4ThreeVector e_perpend(0., 0., 1.);
        if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
        G4ThreeVector e_paralle = e_perpend.cross(kphoton);

        G4ThreeVector polar = std::cos(angle) * e_paralle
                              + std::sin(angle) * e_perpend;
        particle->SetPolarization(polar.x(), polar.y(), polar.z());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


