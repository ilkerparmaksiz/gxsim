//
// Created by ilker on 3/24/23.
//

#ifndef CRAB_OPWLS2_HH
#define CRAB_OPWLS2_HH
#include "G4OpWLS2.hh"
#include "S2Photon.hh"

class OpWLS2:public G4OpWLS2{
public:
   using G4OpWLS2::G4OpWLS2;
   virtual G4bool IsApplicable(const G4ParticleDefinition &aParticleType) override;
};

#endif //CRAB_OPWLS_HH

inline G4bool OpWLS2::IsApplicable(const G4ParticleDefinition &aParticleType) {
    if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
        return true;
    }
    return false;

}