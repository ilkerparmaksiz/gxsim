//
// Created by argon on 2/6/24.
//

#ifndef EXAMPLEB1_PHOTONSD_HH
#define EXAMPLEB1_PHOTONSD_HH
#include "G4VSensitiveDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "PhotonHit.hh"
#include <G4String.hh>
#include <G4Types.hh>
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PhotonSD : public G4VSensitiveDetector{
public:
    PhotonSD(G4String);
    ~PhotonSD() = default;
    void Initialize(G4HCofThisEvent*) final;
    G4bool ProcessHits(G4Step*, G4TouchableHistory* ROhist) final;
    void OpticksHits();
    void EndOfEvent(G4HCofThisEvent* hitCollection) final;
    G4int GetOpticksHits(){ return fOpticksHits;};
    G4int GetGEANT4Hits(){return fGeant4Hits;};

private:
    PhotonHitsCollection* fPhotonHitsCollection{ 0 };
    G4int fHCID{ 0 };
    static constexpr G4double hc = (h_Planck * c_light) / (CLHEP::eV * CLHEP::nm);
    inline G4double etolambda(G4double E)
    {
        // input photon energy in eV
        // return wavelength in nm:
        // lambda = h c/e
        //
        return hc / E;
    }

    G4int fGeant4Hits;
    G4int fOpticksHits;
};


#endif //EXAMPLEB1_PHOTONSD_HH
