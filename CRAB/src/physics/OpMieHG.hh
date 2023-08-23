//
// Created by ilker on 8/24/23.
//

#ifndef CRAB_OpMie_HH
#define CRAB_OpMie_HH
#include "G4OpMieHG.hh"
#include "S2Photon.hh"

class OpMieHG:public G4OpMieHG{
public:
    using G4OpMieHG::G4OpMieHG;
    virtual G4bool IsApplicable(const G4ParticleDefinition &aParticleType) override;
};

#endif //CRAB_OPWLS_HH

inline G4bool OpMieHG::IsApplicable(const G4ParticleDefinition &aParticleType) {
    if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
        return true;
    }
    return false;

}