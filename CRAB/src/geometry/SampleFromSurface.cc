//
// Created by argon on 8/15/23.
//

#include "SampleFromSurface.hh"
#include "G4PVPlacement.hh"
using namespace CLHEP;

SampleFromSurface::SampleFromSurface ():percent_(0.85){
    // Initialize the Map Pointer
    SamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
    TranslatedSamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
}
SampleFromSurface::~SampleFromSurface(){}

// Transforms the sample points
void SampleFromSurface::FaceTransform(const G4VPhysicalVolume* tr) {

  // Fill variables
  const G4String key=tr->GetName();
  G4int TotalPoints= SamplePoints_->at(key).size();

  std::vector<G4ThreeVector> Samples=SamplePoints_->at(key);
  std::vector<G4ThreeVector> TranslatedSamples;
  const G4RotationMatrix ro=*tr->GetRotation();
  const G4ThreeVector translation = tr->GetTranslation();
  G4ThreeVector tVector;

  // Loop Through and apply the transformation
  for (int i=0;i<TotalPoints;i++){
      tVector=translation+(Samples.at(i)*mm);
      //tVector.transform(ro);
      TranslatedSamples.push_back(tVector);
  }

  // Add Translated samples for Generation
    if(TranslatedSamplePoints_->find(tr->GetName())== TranslatedSamplePoints_->end()){
        TranslatedSamplePoints_->insert(std::pair<G4String,std::vector<G4ThreeVector>>(tr->GetName(),TranslatedSamples));
    } else  {
        G4String str;
        str=tr->GetName()+"--> This Key Exists !";
        G4Exception("[SampleFromSurface/FaceTransform]",tr->GetName(),G4ExceptionSeverity::JustWarning,str);
    }



}

void SampleFromSurface::SampleFromFacet(G4String n1, G4TessellatedSolid *solid1) {
    G4int initial;
    G4int TotalFacets=solid1->GetNumberOfFacets();
    if(solid1== nullptr or TotalFacets==0)G4Exception("[SampleFromFacet]","solid",G4ExceptionSeverity::FatalException,"Solid is not Defined or No Facets found");

    if(percent_==0) {
        G4Exception("[SampleFromFacet]","Percent",G4ExceptionSeverity::JustWarning,"Percent is not set! Will Sample from Whole Object");
        initial=0;
    }
    else initial=round(TotalFacets*percent_);

    G4cout <<"Running Facet Sampler For --> "<< n1 <<G4endl;

    std::vector<G4ThreeVector> points;

    // Collect facets into a vector
    for (int i=initial;i<TotalFacets; i++) points.push_back(solid1->GetFacet(i)->GetPointOnFace());

    // Add it to sample
    if(SamplePoints_== nullptr) SamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();


    if(SamplePoints_->find(n1)== SamplePoints_->end()){
        SamplePoints_->insert(std::pair<G4String,std::vector<G4ThreeVector>>(n1,points));
    } else  {
        G4String str;
        str=n1+"--> This Key Exists !";
        G4Exception("[SampleFromFacet]",n1,G4ExceptionSeverity::JustWarning,str);
    }
}
void SampleFromSurface::SampleFromFacet(G4TessellatedSolid *solid1) {
    G4int initial;
    G4int TotalFacets=solid1->GetNumberOfFacets();
    if(solid1== nullptr or TotalFacets==0)G4Exception("[SampleFromFacet]","solid",G4ExceptionSeverity::FatalException,"Solid is not Defined or No Facets found");

    if(percent_==0) {
        G4Exception("[SampleFromFacet]","Percent",G4ExceptionSeverity::JustWarning,"Percent is not set! Will Sample from Whole Object");
        initial=0;
    }
    else initial=round(TotalFacets*percent_);


    G4cout <<"Running Facet Sampler For --> "<< solid1->GetName() <<G4endl;

    std::vector<G4ThreeVector> points;

    // Collect facets into a vector
    for (int i=initial;i<TotalFacets; i++) points.push_back(solid1->GetFacet(i)->GetPointOnFace());

    // Add it to sample
    if(SamplePoints_== nullptr) SamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();


    if(SamplePoints_->find(solid1->GetName())== SamplePoints_->end()){
        SamplePoints_->insert(std::pair<G4String,std::vector<G4ThreeVector>>(solid1->GetName(),points));
    } else  {
        G4String str;
        str=solid1->GetName()+"--> This Key Exists !";
        G4Exception("[SampleFromFacet]",solid1->GetName(),G4ExceptionSeverity::JustWarning,str);
    }
}




