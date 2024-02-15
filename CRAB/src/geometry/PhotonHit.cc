//
// Created by argon on 2/6/24.
//

#include "PhotonHit.hh"
// Geant4 headers
#include "PhotonHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <G4Point3D.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>


template <class Type>
class G4Allocator;
G4ThreadLocal G4Allocator<PhotonHit>* PhotonHitAllocator = nullptr;

PhotonHit::PhotonHit()
        : G4VHit()
{}

PhotonHit::PhotonHit(unsigned iid, unsigned ipid, G4double iwavelength,
                     G4double itime, G4ThreeVector iposition,
                     G4ThreeVector idirection, G4ThreeVector ipolarization)
        : G4VHit()
{
    fid           = iid;
    fpid          = ipid;
    fwavelength   = iwavelength;
    ftime         = itime;
    fposition     = iposition;
    fdirection    = idirection;
    fpolarization = ipolarization;
}

PhotonHit::PhotonHit(const PhotonHit& right)
        : G4VHit()
{
    fid           = right.fid;
    fpid          = right.fpid;
    fwavelength   = right.fwavelength;
    ftime         = right.ftime;
    fposition     = right.fposition;
    fdirection    = right.fdirection;
    fpolarization = right.fpolarization;
}

const PhotonHit& PhotonHit::operator=(const PhotonHit& right)
{
    fid           = right.fid;
    fpid          = right.fpid;
    fwavelength   = right.fwavelength;
    ftime         = right.ftime;
    fposition     = right.fposition;
    fdirection    = right.fdirection;
    fpolarization = right.fpolarization;
    return *this;
}

G4bool PhotonHit::operator==(const PhotonHit& right) const
{
    return (this == &right) ? true : false;
}

void PhotonHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fposition);
        circle.SetScreenSize(2.);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1., 0., 0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}