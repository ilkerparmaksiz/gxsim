#include <iostream>
#include "HeedModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include <stdio.h>
#include "DriftLineTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4EventManager.hh"
#include "G4VVisManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include "G4AutoLock.hh"
namespace{G4Mutex aMutex = G4MUTEX_INITIALIZER;}

const static G4double torr = 1. / 760. * atmosphere;


HeedModel::HeedModel(G4String modelName, G4Region* envelope,DetectorConstruction* dc,GasBoxSD* sd)
  : G4VFastSimulationModel(modelName, envelope), detCon(dc), fGasBoxSD(sd)	{
  fhexrad[1] = 1;
  fhexrad[3] = 7;
  fhexrad[5] = 19; 
  fhexrad[7] = 37;
  fhexrad[9] = 61;
}

HeedModel::~HeedModel() {}

//Method called when a particle is created, checks if the model is applicable for this particle
G4bool HeedModel::IsApplicable(const G4ParticleDefinition& particleType) {
  G4String particleName = particleType.GetParticleName();
  //  G4cout << "HeedModel::IsApplicable(): particle " << particleName  << " allowed? " << FindParticleName(particleName) << G4endl;
  return FindParticleName(particleName);
}

//Method called in every step: checks if the conditions of the particle are met. If true the DoIt-method is called
G4bool HeedModel::ModelTrigger(const G4FastTrack& fastTrack) {
  G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4String particleName =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

  //  G4cout << "HeedModel::ModelTrigger(): particle, ekin " << particleName << ", " << ekin << " allowed? " << FindParticleNameEnergy(particleName, ekin / keV) << G4endl;
  
  return FindParticleNameEnergy(particleName, ekin / keV);
}

//Implementation of the general model, the Run method, calles at the end is specifically implemented for the daughter classes
void HeedModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

  G4ThreeVector dir = fastTrack.GetPrimaryTrack()->GetMomentumDirection();

  G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition();

  G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
  G4String particleName =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

  //  std::cout << "HeedModel::DoIt: track ke [keV], position " << ekin/keV << ",  " << worldPosition[0] << ", " << worldPosition[1] << ", " << worldPosition[2] << std::endl;
  Run(fastStep, fastTrack, particleName, ekin/keV, time, worldPosition.x() / CLHEP::cm,
      worldPosition.y() / CLHEP::cm, worldPosition.z() / CLHEP::cm,
      dir.x(), dir.y(), dir.z());

}

//Checks if the particle is in the list of particle for which the model is applicable (called by IsApplicable)
G4bool HeedModel::FindParticleName(G4String nameloc) {
  MapParticlesEnergy::iterator it;
  it = fMapParticlesEnergy.find(nameloc);
  if (it != fMapParticlesEnergy.end()) {
    return true;
  }
  return false;
}

//Checks if the energy condition of the particle is in the list of conditions for which the model shoould be triggered (called by ModelTrigger)
G4bool HeedModel::FindParticleNameEnergy(G4String name2,
                                             double ekin_keV) {
  MapParticlesEnergy::iterator it;
//  it = fMapParticlesEnergy->find(name2);
  for (it=fMapParticlesEnergy.begin(); it!=fMapParticlesEnergy.end();++it) {
    if(it->first == name2){
      EnergyRange_keV range = it->second;
      if (range.first <= ekin_keV && range.second >= ekin_keV) {
        return true;
      }
    }
  }
  return false;
}

//Initialize the Garfield++ related geometries and physics/tracking mechanisms, this is specific for each use and should be re-implemented entirely
void HeedModel::InitialisePhysics(){
  if(G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::workerRM ||
     G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::sequentialRM ){ // EC added, 9-May-2019.
    makeGas();
      
    buildBox();
    
    BuildCompField();
    
    BuildSensor();
    
    SetTracking();
    
    if(fVisualizeChamber) CreateChamberView();
    if(fVisualizeSignal) CreateSignalView();
    if(fVisualizeField) CreateFieldView();
  }
}

// Gas intialization (see Garfield++ documentation)
void HeedModel::makeGas(){
  fMediumMagboltz = new Garfield::MediumMagboltz();
  double pressure = detCon->GetGasPressure()/torr;
  double temperature = detCon->GetTemperature()/kelvin;
  double arPerc = detCon->GetArgonPercentage();
  double ch4Perc = detCon->GetCH4Percentage();
  std::cout <<  "HeedModel::makeGas: pressure to magboltz is " << pressure << " torr." << std::endl;
  // Another line to comment in/out !   EC, 26-Oct-2021.
  //  fMediumMagboltz->SetComposition("ar", arPerc, "ch4", ch4Perc);
  //fMediumMagboltz->SetComposition("ar", 90, /* "n2", 99, */ "ch4", 10/* "o2", 1*/);
  fMediumMagboltz->SetComposition( "n2", 99,  "o2", 1);
  //  fMediumMagboltz->SetComposition("n2", 100);
  fMediumMagboltz->SetTemperature(temperature);
  fMediumMagboltz->SetPressure(pressure); 
  fMediumMagboltz->EnableDebugging(); 
  fMediumMagboltz->Initialise(true);
  fMediumMagboltz->DisableDebugging();
  
  G4cout << gasFile << G4endl;
  const std::string path = getenv("GARFIELD_HOME");
  G4AutoLock lock(&aMutex);
  if(ionMobFile!="")
    fMediumMagboltz->LoadIonMobility(path + "/Data/" + ionMobFile);
  if(gasFile!=""){
    const bool useLog = true; // to use Logarithmic spacing.
    int nwires = fhexrad.find(detCon->GetNumHexes())->second;
    // After file is successfully created/written, comment out next 3 lines. Also remember to have correct mixture in DetectorConstruction.cc and gasfile name in *.mac.
    // And also note comment above about ar/ch4 vs n2/o2 AND correct gasfile name in .mac file !!

    //    fMediumMagboltz->SetFieldGrid(1.,600000.,20,useLog); // In P10 E/p = 7000000/10 causes Magboltz to crash. EC, 5-Nov-2021.
    /*
    fMediumMagboltz->SetFieldGrid(10.,200000.,20,useLog); // In 99:1 N2:O2 E/p = 3000000/3 causes Magboltz to crash. EC, 5-Nov-2021.`
    fMediumMagboltz->GenerateGasTable(5); // num of 1E7 collisions to do ... something.
    fMediumMagboltz->WriteGasFile(gasFile.c_str());
    */
    fMediumMagboltz->LoadGasFile(gasFile.c_str());


  }
}
  
//Geometry (see Garfield++ documentation)
void HeedModel::buildBox(){
  geo = new Garfield::GeometrySimple();
  double wireR = (80.0E-6*100)/CLHEP::cm;
  // Somehow, this lays everything along the y-axis, and despite AddWire(), etc taking only x,y coordinates it's presumed 
  // that the wires are laid along the y-axis with no explicit rotation required? EC, 11-May-2019. Yes.

  box = new Garfield::SolidTube(0., 0., 0., 0. /*wireR*/,(detCon->GetGasBoxR())/CLHEP::cm,((detCon->GetGasBoxH())*0.5)/CLHEP::cm,0.,1.,0.);

  std::cout << "HeedModel::buildBox() putting cylinder of radius/halfLength " << (detCon->GetGasBoxR())/CLHEP::cm << "/" << (detCon->GetGasBoxH()*0.5)/CLHEP::cm << " [cm] along Garfield y-axis." << std::endl;

  geo->AddSolid(box, fMediumMagboltz);

  std::cout << " .....   now forcing density for the Garfield gas medium -- via P, T." << std::endl;
  geo->GetMedium(0.,0.,0.)->SetPressure(detCon->GetGasPressure()/torr);
  geo->GetMedium(0.,0.,0.)->SetTemperature(detCon->GetTemperature()/kelvin);
  // pull out the density to see if our pressure-setting in .mac file is in fact effected down in Garfield. EC, 28-Oct-2021.
  double ggeodens = geo->GetMedium(0.,0.,0.)->GetMassDensity();
  std::cout << "HeedModel::buildBox(): Garfield mass density [g/cm3] is: " << ggeodens << std::endl;
}

//Construction of the electric field (see Garfield++ documentation)
void HeedModel::BuildCompField(){
    // Switch between IROC and OROC.
    // y-axis gap between rows of wires [cm]
  const double gap = 0.2;
    
    // y coordinates of the wires [cm]

    const double yc = 2. * gap;       // cathode
    const double yg = 2. * gap + 0.3; // gate
    // Periodicity (wire spacing)
    const double period = 0.25;
    const int nRep = 2;
    int nRepC = 10; // randomly chosen number for each side of hex, EC    
    int nRepG = 40; // randomly chosen number for guard ring, EC

    // For 1 Hex I want no cathode/guard wires, just the tube. EC, 14-Oct-2021
    if (detCon->GetNumHexes()==1) {
	nRepC = 0;
	nRepG = 0;
      }

    
    const double dc = period;
    const double dg = period / 2;
    
    // Wire diameters [cm]
    const double dSens = 0.0020;
    const double dCath = 0.0075;
    
    comp = new Garfield::ComponentAnalyticField();
    comp->SetGeometry(geo);
    
    // make 1 hexagon, then repeat it 1, 7, 19, 37, ..., times.

    double nRad = detCon->GetNumHexes(); // read this from config. number of hexes along horizontal diameter (1, 3, 5, 7 ...)

    if (fhexrad.find(nRad) == fhexrad.end()) {
      std::cout << "Throwing, " + std::to_string(nRad) + " not in fhexrad map." << std::endl;
      throw;
    }

    // main hex cathode wires
    double rC = (detCon->GetGasBoxR()/CLHEP::cm - 0.002/CLHEP::cm)*0.95/nRad ; // 95% of distance to 2mm from guard ring wires 
    double xstart = -(nRad-1)*rC ; double ystart = 0.0; 
    double xpt = 0;  double ypt = rC*2./sqrt(3) + rC/2.;  // center point of neighbor hex up, and right
    double ytop = rC*2./sqrt(3); // y coord of top of this hex
    double xhex(0.0); double yhex(0.0);
    int mr(0);
    int inc(0);
    bool newrow(true);
    int rowstop(0);
    for (int nhexes = 0; nhexes<fhexrad.find(nRad)->second; ++nhexes) {

      if (newrow) { // leftmost in row
	xhex = xstart ; yhex = ystart;  newrow = false;
      }
      else if (nhexes<fhexrad.find(nRad)->first) { // on middle row
	xhex += 2*rC ; yhex = 0.0; 
	if (nhexes == fhexrad.find(nRad)->first-1) {
	  mr=1; // finished middle row
	  inc = fhexrad.find(nRad)->first - 1;
	  newrow = true; rowstop = fhexrad.find(nRad)->first + inc;
	  xstart += rC; ystart+= ypt;
	}
      }
      else if (mr%2) { // just finished middle row or a bottom row, now on a top row
	xhex += 2*rC; yhex = ystart;
	if (nhexes == rowstop-1) {
	  mr++; // finished top row
	  newrow = true; rowstop += inc ; // don't decrement rowstop on this top row.
	  ystart-= mr*(ypt);
	}
      }
      else if (mr!=0 && mr%2==0) { // just finished top row, on a bottom row
	xhex += 2*rC; yhex = ystart;
	if (nhexes == rowstop-1) {
	  mr++; // finished bot row
	  newrow = true; inc--; rowstop += inc ; // to force one less hex for each new top row.
	  xstart += rC; ystart+= mr*(ypt);
	}

      }
      // sense wire. 
      // label the middle wire s1 and the last one s2. the rest are generically known as s. 
	std::string sensorName("s");
	if (xhex==0.0 && yhex==0.0) sensorName = "s1";
	if (nhexes == fhexrad.find(nRad)->second-1) sensorName="s2";

	std::cout << "HeedModel::BuildCompField(): Adding anode wire " << sensorName << " of potential " << vAnodeWires << " at " << xhex << ", " << yhex << ". Length [cm] " << detCon->GetGasBoxH()*0.5*4/CLHEP::cm << std::endl;
	comp->AddWire(0.0/CLHEP::cm + xhex, 0.0/CLHEP::cm + yhex, dSens, vAnodeWires, sensorName, detCon->GetGasBoxH()*0.5/CLHEP::cm);
	// bump xs and ys here.

  

	// y axis sides
	double xside = rC + xhex;
	double yside = -rC/2 + yhex; 
	double spacing = rC/nRepC; 

	for (int i = 0; i < nRepC; ++i) {
	  yside+=spacing;
	  comp->AddWire( xside, yside, dCath, vCathodeWires, "c");
	  comp->AddWire(-xside, yside, dCath, vCathodeWires, "c");
	}

	// the upper left side, and its 3 reflection sides
	double xedge = -rC + xhex; // l.l. x coord of top-angle line.
	double yedge = rC/2 + yhex ; 
	// now for the 4 diagonal sides of hexagram
	double spacingx = (xpt - (xedge - xhex))/nRepC; 
	double spacingy = (ytop - (yedge - yhex))/nRepC; 
	for (int i = 0; i < nRepC; ++i) {
	  xedge+=spacingx;
	  yedge+=spacingy;
	  comp->AddWire( xedge, yedge, dCath, vCathodeWires, "c");
	  comp->AddWire(-xedge, yedge, dCath, vCathodeWires, "c"); 
	  comp->AddWire( xedge,-yedge, dCath, vCathodeWires, "c");
	  comp->AddWire(-xedge,-yedge, dCath, vCathodeWires, "c");
	}

    } // end of loop over all hexes

    // guard ring wires. These are the size of cathode wires, and at the sense wire potential.
    double rG = detCon->GetGasBoxR()/CLHEP::cm - 0.002/CLHEP::cm ; // wires offset 2 mm from outer radius of tube
    for (int i = 0; i < nRepG; ++i) {
      double x = rG * cos(2*3.14159*i/nRepG); 
      double y = rG * sin(2*3.14159*i/nRepG); 
      comp->AddWire(x, y, dCath, vAnodeWires, "guard");
    }
    
    double gnd(0.0);
    //    comp->AddTube(detCon->GetGasBoxR()/CLHEP::cm-0.001/CLHEP::cm, gnd, 0, "tubewall");
    comp->AddPlaneY((detCon->GetGasBoxH()*0.5)/CLHEP::cm, gnd, "topcap");
    comp->AddPlaneY(-(detCon->GetGasBoxH()*0.5)/CLHEP::cm, gnd, "bottomcap");

}

//Build sensor (see Garfield++ documentation)
void HeedModel::BuildSensor(){

  // X-width [cm?] of drift simulation will cover between +/- axis_x... rough wire separation, EC.
  const double axis_x = detCon->GetGasBoxR()/detCon->GetNumHexes()/CLHEP::cm;  
  // Y-width of drift simulation will cover between +/- axis_y
  const double axis_z = detCon->GetGasBoxR()/detCon->GetNumHexes()/CLHEP::cm;  
  const double axis_y = detCon->GetGasBoxH()*0.5/CLHEP::cm;

  fSensor = new Garfield::Sensor();
  //  fSensor->EnableDebugging();
  fSensor->AddComponent(comp);
  //  fSensor->SetTimeWindow(0.,2500.,100); //Lowest time [ns], time bins [ns], number of bins
  //  fSensor->SetTimeWindow(0.,5.,100); //Lowest time [ns], time bins [ns], number of bins
  fNumbins = 100;
  fBinsz = 300.; //nsec
  fSensor->SetTimeWindow(0.,fBinsz,fNumbins*fBinsz); //Lowest time [ns], time bins [ns], number of bins

  //comp->AddReadout("s1");
  comp->AddReadout("s2");

  fSensor->SetArea(-axis_x, -axis_y, -axis_z, axis_x, axis_y, axis_z);

  //fSensor->AddElectrode(comp,"s1");
  fSensor->AddElectrode(comp,"s2");

  // Sensor E field check belongs back here. Temporarily moved up to comp.

  std::cout << "HeedModel::BuildSensor() Is geometry a tube?: " << geo->GetSolid(0)->IsTube() << std::endl;  
  std::cout << "HeedModel::BuildSensor() Geometry halflength, radius: " << geo->GetSolid(0)->GetHalfLength() << ", " << geo->GetSolid(0)->GetRadius() << std::endl;
  
}

//Set which tracking mechanism to be used: Runge-kutta, Monte-Carlo or Microscopic (see Garfield++ documentation)
void HeedModel::SetTracking(){
  if(driftRKF){
    fDriftRKF = new Garfield::DriftLineRKF();
    fDriftRKF->SetSensor(fSensor);
    fDriftRKF->EnableDebugging();
  }
  else if(trackMicro){
    fAvalanche = new Garfield::AvalancheMicroscopic();
    fAvalanche->SetSensor(fSensor);
    fAvalanche->EnableSignalCalculation();
  }
  else{  
    fDrift = new Garfield::AvalancheMC();
    fDrift->SetIons();
    fDrift->SetSensor(fSensor);
    fDrift->EnableSignalCalculation();
    fDrift->SetTimeSteps(0.05); // nsec, per example
    fDrift->SetDistanceSteps(2.e-2); // cm, 10x example
    //    fDrift->EnableDebugging(); // way too much information. EC, 2-Nov-2021.
    fDrift->DisableDebugging();
    if(createAval) fDrift->EnableAttachment();
    else fDrift->DisableAttachment();

  }
  fTrackHeed = new Garfield::TrackHeed();
  fTrackHeed->SetSensor(fSensor);
  fTrackHeed->SetParticle("e-");
  fTrackHeed->EnableDeltaElectronTransport();

}

// Set some visualization variables to see tracks and drift lines (see Garfield++ documentation)
void HeedModel::CreateChamberView(){
  char str[30];
  strcpy(str,name);
  strcat(str,"_chamber");
  fChamber = new TCanvas(str, "Chamber View", 700, 700);
  cellView = new Garfield::ViewCell();
  cellView->SetComponent(comp);
  cellView->SetCanvas(fChamber);
  cellView->Plot2d();
  fChamber->Update();
  char str2[30];
  strcpy(str2,name);
  strcat(str2,"_chamber.pdf");
  fChamber->Print(str2);
//  gSystem->ProcessEvents();
  cout << "CreateCellView()" << endl;
  
  viewDrift = new Garfield::ViewDrift();
  viewDrift->SetCanvas(fChamber);
  viewDrift->SetArea(-4.,-4.,+4.,+4.,-10.,+10.);
  if(driftRKF) fDriftRKF->EnablePlotting(viewDrift);
  else if(trackMicro) fAvalanche->EnablePlotting(viewDrift);
  else fDrift->EnablePlotting(viewDrift);
  fTrackHeed->EnablePlotting(viewDrift);

}

//Signal plotting (see Garfield++ documentation)
void HeedModel::CreateSignalView(){
  char str[30];
  strcpy(str,name);
  strcat(str,"_signal");
  fSignal = new TCanvas(str, "Signal on the wire", 700, 700);
  viewSignal = new Garfield::ViewSignal();
  viewSignal->SetCanvas(fSignal);

}

//Electric field plotting (see Garfield++ documentation)
void HeedModel::CreateFieldView(){
  char str[30];
  strcpy(str,name);
  strcat(str,"_efield");
  fField = new TCanvas(name, "Electric field", 700, 700);
  viewField = new Garfield::ViewField();
  viewField->SetCanvas(fField);
  viewField->SetArea(-4.,-4.,+4.,+4.);
  viewField->SetComponent(comp);
  viewField->SetNumberOfContours(40);
  viewField->PlotContour("e");
  fField->Update();
  char str2[30];
  strcpy(str2,name);
  strcat(str2,"_efield.pdf");
  fField->Print(str2);

  int nbinsxz(50); int nbinsy(100); double binedgesxz(4); double binedgesy(10);
  
  TFile* f3d = new TFile("Efield3d.root","recreate");
  TH3F* e3d = new TH3F("E3D","E3Dname",nbinsxz,-binedgesxz,binedgesxz,nbinsy,-binedgesy,binedgesy,nbinsxz,-binedgesxz,binedgesxz);
  TH2F* e2d = new TH2F("E2D","E2Dname",nbinsxz,-binedgesxz,binedgesxz,nbinsy,-binedgesy,binedgesy);

    	    int status(0); Garfield::Medium* medium;
	    double ex,ey,ez;
	    double xmin(+12),xmax(-12),ymin(+12),ymax(-12),zmin(+12),zmax(-12);

	    int nbiny(0);
	    for (double yh=-binedgesy; yh<binedgesy; yh+=2*binedgesy/nbinsy) {
	      nbiny++;
	      int nbinz(0);
	      for (double zr=-binedgesxz; zr<binedgesxz; zr+=2*binedgesxz/nbinsxz) {
		nbinz++;
		int nbinx(0);
		for (double xr=-binedgesxz; xr<binedgesxz; xr+=2*binedgesxz/nbinsxz) {
		  nbinx++;
		  fSensor->ElectricField(xr, yh, zr, ex, ey, ez, medium, status);
		  //		  comp->ElectricField(xr, yh, zr, ex, ey, ez, medium, status);
		  //		  std::cout << "Filling TH3F bin " << xr << "," << yh << "," << zr << std::endl;
		  e3d->SetBinContent(nbinx-1,nbiny-1,nbinz-1,pow(ex*ex+ey*ey+ez*ez,0.5));
		  //		  e2d->SetBinContent(pow((xr+binedgesxz/nbinsxz)*(xr+binedgesxz/nbinsxz)+(zr+binedgesxz/nbinsxz)*(zr+binedgesxz/nbinsxz),0.5),yh+binedgesy/nbinsy,pow(ex*ex+ey*ey+ez*ez,0.5));
		  if (pow(ex*ex+ey*ey+ez*ez,0.5)>0.0001 ) {
		    if (xr<xmin) xmin=xr;
		    if (xr>xmax) xmax=xr;
		    if (zr<zmin) zmin=zr;
		    if (zr>zmax) zmax=zr;
		    if (yh<ymin) ymin=yh;
		    if (yh>ymax) ymax=yh;
	  
		  }
		}}}

	    std::cout << "xmin,xmax; ymin,ymax; zmin,zmax [cm]"<< xmin <<","<<xmax<<"; "<<ymin<< ","<<ymax<<"; "<<zmin<<"," <<zmax  << std::endl;

  
  f3d->Write();
  f3d->Close();
  
}

// Drift the electrons from point of creation towards the electrodes (This is common for both models, i.e. HeedDeltaElectron and HeedModel) (see Garfield++ documentation)
void HeedModel::Drift(double x, double y, double z, double t){
    if(driftElectrons){
        DriftLineTrajectory* dlt = new DriftLineTrajectory();
        G4TrackingManager* fpTrackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
        fpTrackingManager->SetTrajectory(dlt);
        if(driftRKF){
            fDriftRKF->DriftElectron(x,y,z,t);
            unsigned int n = fDriftRKF->GetNumberOfDriftLinePoints();
            double xi,yi,zi,ti;
            for(unsigned int i=0;i<n;i++){
                fDriftRKF->GetDriftLinePoint(i,xi,yi,zi,ti);
                if(G4VVisManager::GetConcreteInstance() && i % 1000 == 0)
                  dlt->AppendStep(G4ThreeVector(xi*CLHEP::cm,yi*CLHEP::cm,zi*CLHEP::cm),ti);
            }
        }
        else if(trackMicro){
            fAvalanche->AvalancheElectron(x,y,z,t,0,0,0,0);
            unsigned int nLines = fAvalanche->GetNumberOfElectronEndpoints();
            for(unsigned int i=0;i<nLines;i++){
                unsigned int n = fAvalanche->GetNumberOfElectronDriftLinePoints(i);
                double xi,yi,zi,ti;
                for(unsigned int j=0;j<n;j++){
                    fAvalanche->GetElectronDriftLinePoint(xi,yi,zi,ti,j,i);
                    if(G4VVisManager::GetConcreteInstance() && i % 1000 == 0)
                      dlt->AppendStep(G4ThreeVector(xi*CLHEP::cm,yi*CLHEP::cm,zi*CLHEP::cm),ti);
                }
            }
        }
        else{
  	    fDrift->DriftElectron(x,y,z,t); // this is where it all happens: drifting, Townshend gain, Signal placement. See my notebook. EC, 2-Nov2021.

            unsigned int n = fDrift->GetNumberOfDriftLinePoints();
            double xi,yi,zi,ti;
	    //	    std::cout << "HeedModeel::Drift(): avalanchetracking, n DLTs is " << n << std::endl;
            for(unsigned int i=0;i<n;i++){
                fDrift->GetDriftLinePoint(i,xi,yi,zi,ti);
                if(G4VVisManager::GetConcreteInstance() ) //&& i % 1000 == 0)
                  dlt->AppendStep(G4ThreeVector(xi*CLHEP::cm,yi*CLHEP::cm,zi*CLHEP::cm),ti);
            }
        }


    }

}

// Plot the track, only called when visualization is turned on by the user
void HeedModel::PlotTrack(){
    if(fVisualizeChamber){
      G4cout << "PlotTrack" << G4endl;
      viewDrift->Plot(true,false);
      fChamber->Update();
      fChamber->Print("PrimaryTrack.pdf");
    }

  //plot the signal here.	
    if(fVisualizeSignal){
      viewSignal->SetSensor(fSensor);
      char str[30];
      strcpy(str,name);
      char str2[30];
      strcpy(str2,name);

      std::cout << "HeedModel::PlotTrack 0"<< std::endl;
      strcat(str,"_signal1.pdf");
      viewSignal->PlotSignal("s1");
      std::cout << "HeedModel::PlotTrack 1"<< std::endl;
      fSignal->Print(str);
      std::cout << "HeedModel::PlotTrack 2"<< std::endl;

      strcat(str2,"_signal2.pdf");
      std::cout << "HeedModel::PlotTrack 3"<< std::endl;
      viewSignal->PlotSignal("s2");
      std::cout << "HeedModel::PlotTrack 4"<< std::endl;
      fSignal->Print(str2);
      std::cout << "HeedModel::PlotTrack 5"<< std::endl;
    }

}
