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
    TFile *FileIn = TFile::Open("macros/crab2MeV_t0.root", "READ");
    // TFile *FileIn = TFile::Open("build/crab2MeV_noref_S2x10.root", "READ");
    


    // CAM --------------------------------------------------------------------
    
    TTree *tCam = (TTree*)FileIn->Get("ntuple/Camera");
    
    Double_t Cam_X, Cam_Y, Cam_Z, Cam_Ev;
    tCam->SetBranchAddress("X",&Cam_X);
    tCam->SetBranchAddress("Y",&Cam_Y);
    tCam->SetBranchAddress("Z",&Cam_Z);
    tCam->SetBranchAddress("Event",&Cam_Ev);

    double SF = 4.1/10; // Converts from mm to cm. Then to the size of the object. M = 4.1

    //create two histograms
    TH2F *hCamXY = new TH2F("hCameraXY","Camera X vs Y; Object Size X [cm]; Object Size Y [cm]",100, -5*SF, 5*SF, 100, -5*SF, 5*SF) ;

    //read all entries and fill the histograms
    Long64_t nCam = tCam->GetEntries();
    for (Long64_t i=0;i<nCam;i++) {
        tCam->GetEntry(i);
        hCamXY->Fill(-Cam_X*SF, -Cam_Y*SF);
    }

    TCanvas *cCam = new TCanvas();
    hCamXY->Draw("colz");

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

    // TCanvas *cLens = new TCanvas();
    // hLensXY->Draw("colz");


    // TCanvas *cLens2 = new TCanvas();
    // hLensXYT->Draw("lego2");


    // S1 ---------------------------------------------------------------------

    TTree *tS1 = (TTree*)FileIn->Get("ntuple/S1");
    
    Double_t S1_X, S1_Y, S1_Z, S1_Ev;
    tS1->SetBranchAddress("X",&S1_X);
    tS1->SetBranchAddress("Y",&S1_Y);
    tS1->SetBranchAddress("Z",&S1_Z);
    tS1->SetBranchAddress("Event",&S1_Ev);

    //create two histograms
    TH3F *hS1XYZ = new TH3F("hS1XYZ","S1 XYZ; X [mm]; Y [mm]; Z [mm]",40, -36, 36, 40, -36, 36, 40, -220, 220) ;
    TH2F *hS1XY = new TH2F("hS1XY","S1 XY; X [cm]; Y [cm]",80, -3.6, 3.6, 80, -3.6, 3.6) ;

    //read all entries and fill the histograms
    Long64_t nS1 = tS1->GetEntries();
    for (Long64_t i=0;i<nS1;i++) {
        tS1->GetEntry(i);
        hS1XYZ->Fill(S1_X, S1_Y, S1_Z);
        hS1XY->Fill(S1_X/10, S1_Y/10);
    }

    // TCanvas *cS1 = new TCanvas();
    // hS1XYZ->Draw("lego2");
    
    TCanvas *cS12 = new TCanvas();
    hS1XY->Draw("colz");


    // S2 ---------------------------------------------------------------------

    TTree *tS2 = (TTree*)FileIn->Get("ntuple/S2");
    
    Double_t S2_X, S2_Y, S2_Z, S2_t,S2_Ev;
    tS2->SetBranchAddress("X",&S2_X);
    tS2->SetBranchAddress("Y",&S2_Y);
    tS2->SetBranchAddress("Z",&S2_Z);
    tS2->SetBranchAddress("Time",&S2_t);
    tS2->SetBranchAddress("Event",&S2_Ev);

    //create two histograms
    TH1F *hS2T = new TH1F("hS2T","S2 Timing; Time [us]; Entries",500, 40, 120) ;

    //read all entries and fill the histograms
    Long64_t nS2 = tS2->GetEntries();
    
    for (Long64_t i=0;i<nS2;i++) {
        tS2->GetEntry(i);

        hS2T->Fill(S2_t/1000.);

    }

    TCanvas *cS2 = new TCanvas();
    hS2T->Draw("hist");

    // PMT --------------------------------------------------------------------
    
    TTree *tPMT= (TTree*)FileIn->Get("ntuple/PMT");
    
    Double_t PMT_X, PMT_Y, PMT_Z, PMT_Ev;
    tPMT->SetBranchAddress("X",&PMT_X);
    tPMT->SetBranchAddress("Y",&PMT_Y);
    tPMT->SetBranchAddress("Z",&PMT_Z);
    tPMT->SetBranchAddress("Event",&PMT_Ev);

    //create two histograms
    TH2F *hPMTXY = new TH2F("hPMTXY","PMT X vs Y",50, -12.7, 12.7, 50, -12.7, 12.7) ;

    //read all entries and fill the histograms
    Long64_t nPMT= tPMT->GetEntries();
    for (Long64_t i=0;i<nPMT;i++) {
        tPMT->GetEntry(i);
        hPMTXY->Fill(-PMT_X, -PMT_Y);
    }

    TCanvas *cPMT= new TCanvas();
    hPMTXY->Draw("colz");





}
