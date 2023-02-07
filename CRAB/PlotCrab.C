// A root macro to plot the information stored in the CRAB output files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TGeoTube.h>

void PlotCrab(){



    // TFile *FileIn = TFile::Open("build/crab50MeVeminus_noref.root", "READ");
    // TFile *FileIn = TFile::Open("build/crab50MeVeminus_withref.root", "READ");
    TFile *FileIn = TFile::Open("build/crab2MeV.root", "READ");
    // TFile *FileIn = TFile::Open("build/crab2MeV_noref_S2x10.root", "READ");


    // CAM --------------------------------------------------------------------
    // Code uses PMT, so sticking with that convention. Really we mean camera.

    TTree *tPMT = (TTree*)FileIn->Get("ntuple/PMT");
    
    Double_t PMT_X, PMT_Y, PMT_Z, PMT_Ev;
    tPMT->SetBranchAddress("X",&PMT_X);
    tPMT->SetBranchAddress("Y",&PMT_Y);
    tPMT->SetBranchAddress("Z",&PMT_Z);
    tPMT->SetBranchAddress("Event",&PMT_Ev);

    //create two histograms
    TH2F *hPMTXY = new TH2F("hPMTXY","Camera X vs Y",50, -15, 15, 50, -15, 15) ;

    //read all entries and fill the histograms
    Long64_t nPMT = tPMT->GetEntries();
    for (Long64_t i=0;i<nPMT;i++) {
        tPMT->GetEntry(i);
        hPMTXY->Fill(-PMT_X, -PMT_Y);
    }

    TCanvas *cPMT = new TCanvas();
    hPMTXY->Draw("colz");

    // LENS -------------------------------------------------------------------

    TTree *tLens = (TTree*)FileIn->Get("ntuple/Lens");
    
    Double_t Lens_X, Lens_Y, Lens_Z, Lens_Ev, Lens_T;
    tLens->SetBranchAddress("X",&Lens_X);
    tLens->SetBranchAddress("Y",&Lens_Y);
    tLens->SetBranchAddress("Z",&Lens_Z);
    tLens->SetBranchAddress("Event",&Lens_Ev);
    tLens->SetBranchAddress("Time",&Lens_T);

    //create two histograms
    TH2F *hLensXY = new TH2F("hLensXY","Lens X vs Y",10, -5, 5, 10, -5, 5) ;

    // XZ vs T
    TH3F *hLensXYT = new TH3F("hLensXXT","Lens XY vs T",50, 80e3, 130e3, 10, -5, 5, 10, -5, 5 ) ;

    //read all entries and fill the histograms
    Long64_t nLens = tLens->GetEntries();
    for (Long64_t i=0;i<nLens;i++) {
        tLens->GetEntry(i);
        hLensXY->Fill(Lens_X, Lens_Y);
        hLensXYT->Fill(Lens_T, Lens_X, Lens_Y);
    }

    TCanvas *cLens = new TCanvas();
    hLensXY->Draw("colz");


    TCanvas *cLens2 = new TCanvas();
    hLensXYT->Draw("lego2");


    // S1 ---------------------------------------------------------------------

    TTree *tS1 = (TTree*)FileIn->Get("ntuple/S1");
    
    Double_t S1_X, S1_Y, S1_Z, S1_Ev;
    tS1->SetBranchAddress("X",&S1_X);
    tS1->SetBranchAddress("Y",&S1_Y);
    tS1->SetBranchAddress("Z",&S1_Z);
    tS1->SetBranchAddress("Event",&S1_Ev);

    //create two histograms
    TH3F *hS1XYZ = new TH3F("hS1XYZ","S1 XYZ; X [mm]; Y [mm]; Z [mm]",40, -42.5, 42.5, 40, -42.5, 42.5, 40, -220, 220) ;
    TH2F *hS1XY = new TH2F("hS1XY","S1 XY; X [mm]; Y [mm]",40, -42.5, 42.5, 40, -42.5, 42.5) ;

    //read all entries and fill the histograms
    Long64_t nS1 = tS1->GetEntries();
    for (Long64_t i=0;i<nS1;i++) {
        tS1->GetEntry(i);
        hS1XYZ->Fill(S1_X, S1_Y, S1_Z);
        hS1XY->Fill(S1_X, S1_Y);
    }

    TCanvas *cS1 = new TCanvas();
    hS1XYZ->Draw("lego2");
    
    TCanvas *cS12 = new TCanvas();
    hS1XY->Draw("colz");

    


}