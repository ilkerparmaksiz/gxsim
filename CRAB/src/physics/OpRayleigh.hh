//
// Created by ilker on 8/24/23.
//

#ifndef CRAB_OpRayleigh_HH
#define CRAB_OpRayleigh_HH
#include "G4OpRayleigh.hh"
#include "S2Photon.hh"

class OpRayleigh:public G4OpRayleigh{
public:
    using G4OpRayleigh::G4OpRayleigh;
    virtual G4bool IsApplicable(const G4ParticleDefinition &aParticleType) override;
};

#endif //CRAB_OPWLS_HH

inline G4bool OpRayleigh::IsApplicable(const G4ParticleDefinition &aParticleType) {
    if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
        return true;
    }
    return false;

}