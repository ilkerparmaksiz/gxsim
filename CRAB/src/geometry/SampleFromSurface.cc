//
// Created by Ilker Parmaksiz on 8/15/23.
//

#include "SampleFromSurface.hh"
#include "G4PVPlacement.hh"
#include "FileHandling.hh"
#include "Randomize.hh"
#include <random>

using namespace CLHEP;
using namespace util;

 SampleFromSurface::SampleFromSurface (G4String name):frac_(0.135),OverRide(false){
     name_=name;
     SamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
     TranslatedSamplePoints_=new std::map<G4String,std::vector<G4ThreeVector>>();
     file=new filehandler::FileHandling();

     CRABPATH=getenv("CRABPATH");
     if(CRABPATH=="") G4Exception("SampleFromSurface","CRABPATH",FatalException,"CRAB Path is not defined");
     AllFilePath=CRABPATH+"data/"+name_+".txt";
}

SampleFromSurface::~SampleFromSurface(){}

// Transforms the sample points so they will allign with the geometry
void SampleFromSurface::FaceTransform(const G4VPhysicalVolume* tr,const G4VPhysicalVolume * Mother) {
   const G4String key=tr->GetName();
   G4String SingleFilePath=CRABPATH+"data/"+key+".txt";
   // If the File Exist do not worry of transforming it
   if(file->FileCheck(SingleFilePath) and OverRide==false){
       G4String ss;
       ss=key+".txt is exist !";
       G4Exception("[SampleFromSurface]","FaceTransform",JustWarning,ss);
       return;
   }

  // Fill variables
  G4int TotalPoints= SamplePoints_->at(key).size();

  std::vector<G4ThreeVector> Samples=SamplePoints_->at(key);
  std::vector<G4ThreeVector> TranslatedSamples;
  const G4RotationMatrix ro=tr->GetObjectRotationValue();
  G4ThreeVector Mothertranslation = Mother->GetObjectTranslation();
  const G4RotationMatrix MotherRotation = Mother->GetObjectRotationValue();
  G4ThreeVector Daughtertranslation = tr->GetObjectTranslation();
  G4ThreeVector tVector;

  // Loop Through and apply the transformation
  for (int i=0;i<TotalPoints;i++){
      tVector=Daughtertranslation+(Samples.at(i)*mm).transform(ro).transform(MotherRotation)+Mothertranslation;
      TranslatedSamples.push_back(tVector);
  }

  // Add Translated samples for Generation
    if(TranslatedSamplePoints_->find(tr->GetName())== TranslatedSamplePoints_->end()){
        TranslatedSamplePoints_->insert(std::pair<G4String,std::vector<G4ThreeVector>>(tr->GetName(),TranslatedSamples));
    } else {
        G4String str;
        str=tr->GetName()+"--> This Key Exists !";
        G4Exception("[SampleFromSurface/FaceTransform]",tr->GetName(),G4ExceptionSeverity::JustWarning,str);
    }
    // Create a Single File
    //file->SaveToTextFile(SingleFilePath,"X,Y,Z",',',TranslatedSamples);
}
// Sampling from vertex points of the triangles
void SampleFromSurface::SampleFromFacet(G4String n1, G4TessellatedSolid *solid1) {
    G4int initial;
    G4int TotalFacets=solid1->GetNumberOfFacets();
    if(file->FileCheck(AllFilePath) and OverRide==false){
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
    Hep3Vector tmpvect;
    // Collect facets into a vector
    for (int i=initial;i<TotalFacets; i++) {
        //Round up
        tmpvect[0]=round(solid1->GetFacet(i)->GetPointOnFace()[0]*100)/100;
        tmpvect[1]=round(solid1->GetFacet(i)->GetPointOnFace()[1]*100)/100;
        tmpvect[2]=round(solid1->GetFacet(i)->GetPointOnFace()[2]*100)/100;
        points.push_back(tmpvect);
    }

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
    G4String SingleFilePath=CRABPATH+"data/"+solid1->GetName()+".txt";
    if(file->FileCheck(SingleFilePath)  and OverRide==false){
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
    Hep3Vector tmpvect;
    // Collect facets into a vector
    for (int i=initial;i<TotalFacets; i++){
        tmpvect[0]=round(solid1->GetFacet(i)->GetPointOnFace()[0]*100)/100;
        tmpvect[1]=round(solid1->GetFacet(i)->GetPointOnFace()[1]*100)/100;
        tmpvect[2]=round(solid1->GetFacet(i)->GetPointOnFace()[2]*100)/100;
        points.push_back(tmpvect);
    }

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
void SampleFromSurface::SaveAllPointsToOneFile() {
    if(file->FileCheck(AllFilePath) and OverRide==false){
        G4String ss;
        ss=name_+".txt is exist !";
        G4Exception("[SampleFromSurface]","SaveAllPointsToOneFile",JustWarning,ss);
        return;
    }
    if(SamplePoints_->empty()) G4Exception("[SampleFromSurface::SaveAllPointsToOneFile]","SamplePoints",FatalException,"This variable empty !!");

    if(TranslatedSamplePoints_->empty())
        //Assuming no need to call FaceTransform to take account of shifting or rotation
        TranslatedSamplePoints_=SamplePoints_;
    // Create SingleFiles
    // Randomly Distribute Points
    std::random_device frandomDevice;
    std::mt19937 fGenerator(frandomDevice());
    std::vector<G4ThreeVector> All;
     // Break down the map to its key (G4String) and points (std::vector<std::vector<Hep3Vector>>)
     G4String SingleFilePath;
     for (auto &mapValue : *TranslatedSamplePoints_){
           // Save everything to one array
         SingleFilePath=CRABPATH+"data/"+mapValue.first+".txt";
         // Lets Shuffle the Points
         std::shuffle(mapValue.second.begin(),mapValue.second.end(),fGenerator);
         file->SaveToTextFile(SingleFilePath,"X,Y,Z",',',mapValue.second);
         for(int i=0;i<mapValue.second.size();i++) All.push_back(mapValue.second.at(i));
     }

    std::shuffle(All.begin(),All.end(),fGenerator);
    std::shuffle(All.begin(),All.end(),fGenerator);

    // SaveThemTo File
    file->SaveToTextFile(AllFilePath,"X,Y,Z",',',All);
}




