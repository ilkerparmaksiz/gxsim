// A root macro to calculate the overall yields from alpha events
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TGeoTube.h>


void GetYield(TFile* file, std::string TreeName, std::vector<double> &Yield){

    TTree *Tree = (TTree*)file->Get(TreeName.c_str());


    std::vector<double> LY; // Light yield at the camera
    
    Double_t X, Y, Z, Ev;
    Tree->SetBranchAddress("X",&X);
    Tree->SetBranchAddress("Y",&Y);
    Tree->SetBranchAddress("Z",&Z);
    Tree->SetBranchAddress("Event",&Ev);


    double temp_LY = 0;
    double temp_event = -99;

    //read all entries and fill the histograms
    Long64_t N = Tree->GetEntries();
    for (Long64_t i=0; i < N; i++) {
        
        Tree->GetEntry(i);

        if (temp_event == -99){
            temp_event = Ev;
            temp_LY +=1;
            continue;
        }

        if (temp_event != Ev){
            temp_event = Ev;
            Yield.push_back(temp_LY*0.18);
            std::cout << temp_LY << std::endl;
            temp_LY = 0;
        }
        else {
             temp_LY+=1;
        }
    }

}


void CalcYields(){

    TFile *FileIn = TFile::Open("macros/alpha_merge.root", "READ");


    // CAM --------------------------------------------------------------------
    std::vector<double> camYields;
    GetYield(FileIn, "ntuple/Camera", camYields);

    TH1D *hLY = new TH1D("hLY_Camera", "Camera LY; LY; Entries", 200, 0 , 5000);

    for (unsigned int i = 0; i < camYields.size(); i++){
        hLY->Fill(camYields.at(i));
    }

    TCanvas *cCamLY = new TCanvas();
    hLY->Draw("hist");


    // PMT --------------------------------------------------------------------
    std::vector<double> PMTYields;
    GetYield(FileIn, "ntuple/PMT", PMTYields);

    TH1D *hLY_PMT = new TH1D("hLY_PMT", "PMT LY; LY; Entries", 200, 0 , 12000);

    for (unsigned int i = 0; i < PMTYields.size(); i++){
        hLY_PMT->Fill(PMTYields.at(i));
    }

    TCanvas *cPMTLY = new TCanvas();
    hLY_PMT->Draw("hist");

}