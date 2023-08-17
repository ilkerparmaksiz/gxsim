//
// Created by Ilker Parmaksiz on 8/15/23.
//

#include "SampleFromSurface.hh"
#include "G4PVPlacement.hh"
#include "FileHandling.hh"
using namespace CLHEP;
using namespace util;

 SampleFromSurface::SampleFromSurface (G4String name):frac_(0.135),OverRide(true){
     name_=name;
     SamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
     TranslatedSamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
     file=new filehandler::FileHandling();

     std::string crabpath= getenv("CRABPATH");
     if(crabpath=="") G4Exception("SampleFromSurface","CRABPATH",FatalException,"CRAB Path is not defined");
     FilePath=crabpath+"data/"+name_+".txt";
}

SampleFromSurface::~SampleFromSurface(){}

// Transforms the sample points so they will allign with the geometry
void SampleFromSurface::FaceTransform(const G4VPhysicalVolume* tr) {

   // If the File Exist do not worry of transforming it
   if(file->FileCheck(FilePath) and OverRide==false){
       G4String ss;
       ss=name_+".txt is exist !";
       G4Exception("[SampleFromSurface]","FaceTransform",JustWarning,ss);
       return;
   }

  // Fill variables
  const G4String key=tr->GetName();
  G4int TotalPoints= SamplePoints_->at(key).size();

  std::vector<G4ThreeVector> Samples=SamplePoints_->at(key);
  std::vector<G4ThreeVector> TranslatedSamples;
  const G4RotationMatrix ro=tr->GetObjectRotationValue();
  G4ThreeVector translation = tr->GetObjectTranslation();
  G4ThreeVector tVector;

  // Loop Through and apply the transformation
  for (int i=0;i<TotalPoints;i++){
      tVector=translation+(Samples.at(i)*mm).transform(ro);
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

    // Save it to the file
    file->SaveToTextFile(FilePath,"X,Y,Z",',',TranslatedSamples);
}
// Sampling from vertex points of the triangles
void SampleFromSurface::SampleFromFacet(G4String n1, G4TessellatedSolid *solid1) {
    G4int initial;
    G4int TotalFacets=solid1->GetNumberOfFacets();
    if(file->FileCheck(FilePath) and OverRide==false){
        G4String ss;
        ss=name_+".txt is exist !";
        G4Exception("[SampleFromSurface]","SampleFromFacet",JustWarning,ss);
        return;
    }
    if(solid1== nullptr or TotalFacets==0)G4Exception("[SampleFromFacet]","solid",G4ExceptionSeverity::FatalException,"Solid is not Defined or No Facets found");

    if(frac_==0) {
        G4Exception("[SampleFromFacet]","Fraction",G4ExceptionSeverity::JustWarning,"Percent is not set! Will Sample from Whole Object");
        initial=0;
    }
    else initial=round(TotalFacets*(1-frac_));

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
    if(file->FileCheck(FilePath)  and OverRide==false){
        G4String ss;
        ss=name_+".txt is exist !";
        G4Exception("[SampleFromSurface]","SampleFromFacet",JustWarning,ss);
        return;
    }
    if(solid1== nullptr or TotalFacets==0)G4Exception("[SampleFromFacet]","solid",G4ExceptionSeverity::FatalException,"Solid is not Defined or No Facets found");

    if(frac_==0) {
        G4Exception("[SampleFromFacet]","Percent",G4ExceptionSeverity::JustWarning,"Percent is not set! Will Sample from Whole Object");
        initial=0;
    }
    else initial=round(TotalFacets*(1-frac_));

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




