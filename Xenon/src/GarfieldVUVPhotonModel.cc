#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "GarfieldVUVPhotonModel.hh"
#include "G4Region.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"
#include <fstream>
#include "G4TransportationManager.hh"
#include "G4DynamicParticle.hh"
#include "G4RandomDirection.hh"
#include "globals.hh"
#include "MediumMagboltz.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "Medium.hh"
#include "SolidTube.hh"
#include "G4OpticalPhoton.hh"
#include "GarfieldExcitationHit.hh"
#include "GasModelParameters.hh"
#include "DetectorConstruction.hh"
#include "GasBoxSD.hh"
#include "G4ProcessManager.hh"


const static G4double torr = 1. / 760. * atmosphere;

GarfieldVUVPhotonModel::GarfieldVUVPhotonModel(GasModelParameters* gmp, G4String modelName,G4Region* envelope,DetectorConstruction* dc,GasBoxSD* sd) :
		G4VFastSimulationModel(modelName, envelope),detCon(dc),fGasBoxSD(sd) {
	thermalE=gmp->GetThermalEnergy();
	InitialisePhysics();
}

G4bool GarfieldVUVPhotonModel::IsApplicable(const G4ParticleDefinition& particleType) {

  if (particleType.GetParticleName()=="e-") {
      G4cout<<"VUVModel::IsApplic(): True"<<G4endl;
      return true;
  }
      G4cout<<"VUVModel::IsApplic(): False"<<G4endl;
  	return false;
		
		
}

G4bool GarfieldVUVPhotonModel::ModelTrigger(const G4FastTrack& fastTrack){
  G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4cout<<"VUVModel::ModelTrigg(): ekin, ThermalE: " << ekin<< ", "<< thermalE <<G4endl;
  if (ekin<thermalE)
		return true;
  return false;

} 
	
void GarfieldVUVPhotonModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) 
{

  /* 
     This takes all the tracks from the  ionization/conversion electrons created by Degrad, then tracked by G4 from Degrad's simulation 
     of the primary gamma's photoelectric effect in Xe. 
     Each such electron is drifted in the E-field and avalanched, as appropriate. That creates excited Xe atoms. We put one ELM photon
     per excitation of 172/7.2 nm/eV onto the optical photon stack. Geant4 will track those in the normal way. 
     (I think Garfield creates these photons too, but we're throwing 'em out.)

     Note, the weirdness of userHandle, which seems to be called at each excitation by the AvalancheMicroscopic model and fills our
     GarfieldExcitationHitCollection. And which somehow we're allowed to grab here.

     EC, 2-Dec-2021.

   */
  
    G4cout<<"HELLO Garfield 0"<<G4endl;
    ////The details of the Garfield model are implemented here
     fastStep.KillPrimaryTrack();//KILL DEGRAD TRACKS
     G4cout<<"HELLO Garfield 1"<<G4endl;
     garfPos =fastTrack.GetPrimaryTrack()->GetVertexPosition();
         G4cout<<"HELLO Garfield 2"<<G4endl;
     garfTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
     //G4cout<<"GLOBAL TIME "<<G4BestUnit(garfTime,"Time")<<" POSITION "<<G4BestUnit(garfPos,"Length")<<G4endl;


    GenerateVUVPhotons(fastTrack,fastStep,garfPos,garfTime);
    G4cout<<"HELLO Garfield 3"<<G4endl;

}

GarfieldExcitationHitsCollection *garfExcHitsCol;

void GarfieldVUVPhotonModel::GenerateVUVPhotons(const G4FastTrack& fastTrack, G4FastStep& fastStep,G4ThreeVector garfPos,G4double garfTime)
{


	G4double x0=garfPos.getX()*0.1;//Garfield length units are in cm
	G4double y0=garfPos.getY()*0.1;
	G4double z0=garfPos.getZ()*0.1;
	G4double t0=garfTime;
	G4double e0=7.;// starting energy [eV]->I have chose 7 eV because is the energy cut in Degrad
	garfExcHitsCol = new GarfieldExcitationHitsCollection();
	fAvalanche->AvalancheElectron(x0,y0,z0,t0, e0, 0., 0., 0.);

	unsigned int nElastic, nIonising, nAttachment, nInelastic, nExcitation, nSuperelastic;
	fMediumMagboltz->GetNumberOfElectronCollisions(nElastic, nIonising, nAttachment, nInelastic, nExcitation, nSuperelastic);
	
	G4cout<<"NExcitation "<<nExcitation<<G4endl; // This quantity seems to be cumulative over (at least) the event. ... EC, 2-Dec-2021.

	G4int colHitsEntries=garfExcHitsCol->entries();
	G4cout<<"GarfExcHits entries "<<colHitsEntries<<G4endl; // This one is not cumulative.
	
	for (G4int i=0;i<colHitsEntries;i++){
	  GarfieldExcitationHit* newExcHit=new GarfieldExcitationHit();
	  newExcHit->SetPos((*garfExcHitsCol)[i]->GetPos());
	  newExcHit->SetTime((*garfExcHitsCol)[i]->GetTime());
	  fGasBoxSD->InsertGarfieldExcitationHit(newExcHit);
	  // fastStep.SetNumberOfSecondaryTracks(1);	//1 photon per excitation .... Must be commented when I comment below condition too.
	  if(i % (colHitsEntries/10) == 0){ // Need to uncomment this condition, along with one in degradmodel.cc. EC, 2-Dec-2021.
	    G4DynamicParticle VUVphoton(G4OpticalPhoton::OpticalPhotonDefinition(),G4RandomDirection(), 7.2*eV);
	    // Create photons track
	    G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, (*garfExcHitsCol)[i]->GetPos(),(*garfExcHitsCol)[i]->GetTime(),false);
	    G4ProcessManager* pm= newTrack->GetDefinition()->GetProcessManager();
	    //	G4ProcessVectorfAtRestDoItVector = pm->GetAtRestProcessVector(typeDoIt);
	    //G4ProcessVectorfVector gpv = pm->GetAtRestProcessVector(typeGPIL);
	    }						
	}
	delete garfExcHitsCol;
}
// Selection of Xenon exitations and ionizations

void GarfieldVUVPhotonModel::InitialisePhysics(){
	fMediumMagboltz = new Garfield::MediumMagboltz();
	double pressure = detCon->GetGasPressure()/torr;
	double temperature = detCon->GetTemperature()/kelvin;

	fMediumMagboltz->SetTemperature(temperature);
	fMediumMagboltz->SetPressure(pressure);
	fMediumMagboltz->SetComposition("Xe", 100.);

	Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
	// Make a box
	G4double detectorRadius=detCon->GetGasBoxR();//cm
	G4double detectorHalfZ=detCon->GetGasBoxH()*0.5;//cm

	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 0" << std::endl;
	//	Garfield::SolidTube* tube = new Garfield::SolidTube(0.0, detectorHalfZ/CLHEP::cm,0.,0.0, detectorRadius/CLHEP::cm,detectorHalfZ/CLHEP::cm,0.,1.,0.);//Tube oriented in Y'axis (0.,1.,0.,)
	Garfield::SolidTube* tube = new Garfield::SolidTube(0.0, detectorHalfZ/CLHEP::cm*2.0+0.5*detectorHalfZ/CLHEP::cm*2.*0.05,0.,0.0, detectorRadius/CLHEP::cm,detectorHalfZ/CLHEP::cm*0.05,0.,1.,0.);//Tube oriented in Y'axis (0.,1.,0.,)
	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 1" << std::endl;
	// Add the solid to the geometry, together with the medium inside
	geo->AddSolid(tube, fMediumMagboltz);
	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 2" << std::endl;
	// Make a component with analytic electric field
	Garfield::ComponentConstant* componentConstant = new Garfield::ComponentConstant();
	componentConstant->SetGeometry(geo);
	//SetElectricField(const double ex, const double ey, const double ez);

	componentConstant->SetElectricField(0.,-3000.0,0.);
	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 3" << std::endl;

	// Make a sensor
	Garfield::Sensor* sensor = new Garfield::Sensor();
	sensor->AddComponent(componentConstant);

	fAvalanche = new Garfield::AvalancheMicroscopic();
	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 4" << std::endl;

	fAvalanche->SetUserHandleInelastic(userHandle);
		
	fAvalanche->SetSensor(sensor);
	std::cout << "GarfieldVUVPhotonModel::InitPhys(): 5" << std::endl;
  
}

// Selection of Xenon excitations and ionizations
void userHandle(double x, double y, double z, double t, int type, int level,Garfield::Medium * m)
{
	G4ThreeVector Pos;

	if (level > 2 && level < 53){//XENON LEVELS
	
	GarfieldExcitationHit* newExcHit=new GarfieldExcitationHit();
	Pos.setX(x*10);//back to cm to GEANT4
	Pos.setY(y*10);//back to cm to GEANT4
	Pos.setZ(z*10);//back to cm to GEANT4
	newExcHit->SetPos(Pos);
	newExcHit->SetTime(t);
	garfExcHitsCol->insert(newExcHit);
	//If choose to draw change the visualizer from OGL to HepRep in vis.mac file
	//newExcHit->Draw();	
	}// if level
	
	
}// end userhandle	
