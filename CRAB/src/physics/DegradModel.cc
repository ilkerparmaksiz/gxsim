#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "DegradModel.hh"
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
#include "GasModelParameters.hh"
#include "GasBoxSD.hh"
#include "XenonHit.hh"
#include "G4VProcess.hh"
#include "DetectorConstruction.hh"
#include "G4/NESTProc.hh"

DegradModel::DegradModel(GasModelParameters* gmp, G4String modelName, G4Region* envelope,DetectorConstruction* dc, GasBoxSD* sd)
    : G4VFastSimulationModel(modelName, envelope),detCon(dc), fGasBoxSD(sd){
    
    thermalE=gmp->GetThermalEnergy();
    processOccured = false;

    // Get the path to crab
    crab_path = std::getenv("CRABPATH");
    if (crab_path == nullptr) {
        G4Exception("[DegradModel]", "DegradModel()", FatalException,
                    "Environment variable CRABPATH not defined!");
    }

    G4String path(crab_path);

}

DegradModel::~DegradModel() {}

G4bool DegradModel::IsApplicable(const G4ParticleDefinition& particleType) {
    if (particleType.GetParticleName()=="e-")
        return true;
    
    return false;
}

G4bool DegradModel::ModelTrigger(const G4FastTrack& fastTrack) {
    G4int id = fastTrack.GetPrimaryTrack()->GetParentID();
    
    if (id == 1){
        //  also require that only photoelectric effect electrons are tracked here.
        if ( (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("phot") != std::string::npos) ||
             (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("comp") != std::string::npos))
            return true;
    }
  
    return false;

}

void DegradModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

    // Here we start by killing G4's naive little one photo-electron.
    // Then we run degrad with a photon of desired energy. Then we read up all the electrons it produces.
    fastStep.KillPrimaryTrack();
    if(!processOccured){
        G4ThreeVector degradPos =fastTrack.GetPrimaryTrack()->GetVertexPosition();
        G4double degradTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
        G4int KE = int(fPrimPhotonKE/eV);
        const static G4double torr = 1. / 750.062 * bar;
        G4int Press = int(detCon->GetGasPressure()/torr); 
        fastStep.SetPrimaryTrackPathLength(0.0);
        G4cout<<"GLOBAL TIME "<<G4BestUnit(degradTime,"Time")<<" POSITION "<<G4BestUnit(degradPos,"Length")<<G4endl;

        G4int stdout;
        G4int SEED=54217137*G4UniformRand();
        G4String seed = G4UIcommand::ConvertToString(SEED);
        // Note the exact precision in below arguments. The integers gammaKE,xenonP in particular need a ".0" tacked on.
        // G4String degradString="printf \"1,1,3,-1,"+seed+",30000.0,7.0,0.0\n7,0,0,0,0,0\n100.0,0.0,0.0,0.0,0.0,0.0,20.0,900.0\n500.0,0.0,0.0,1,0\n100.0,0.5,1,1,1,1,1,1,1\n0,0,0,0,0,0\" > conditions_Degrad.txt";
        G4String gammaKE(","+std::to_string(KE));
        G4String xenonP(","+std::to_string(Press));
        G4String degradString="printf \"1,1,3,-1,"+seed+gammaKE+".0,7.0,0.0\n7,0,0,0,0,0\n100.0,0.0,0.0,0.0,0.0,0.0,20.0"+xenonP+".0\n500.0,0.0,0.0,1,0\n100.0,0.5,1,1,1,1,1,1,1\n0,0,0,0,0,0\" > conditions_Degrad.txt";
        G4cout << degradString << G4endl;
        stdout=system(degradString.data());
        G4cout << degradString << G4endl;
        const std::string degradpath = std::getenv("DEGRAD_HOME");
        G4cout << degradpath << G4endl;
        std::string exec = "/Degrad < conditions_Degrad.txt";
        std::string full_path = degradpath + exec;
        const char *mychar = full_path.c_str();
        G4cout << mychar << G4endl;
        stdout=system(mychar);
        stdout=system("/convertDegradFile.py");

        GetElectronsFromDegrad(fastStep,degradPos,degradTime);
        // We call Degrad only once, which now that we have the x,y,z location of our primary Xray interaction, re-simulates that interaction. 
        // Note the 5900 in the Degrad config file. The above system() line forces the single Degrad execution. EC, 2-Dec-2021.
        processOccured=true; 
    }

}

void DegradModel::GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime)
{
    G4int eventNumber,Nep, nline, i, electronNumber; //Nep é o numero de e primarios que corresponde ao que o biagi chama de "ELECTRON CLUSTER SIZE (NCLUS)"
    G4double posX,posY,posZ,time,n;
    G4double  posXDegrad,posYDegrad,posZDegrad,timeDegrad;
    G4double  posXInitial=degradPos.getX();
    G4double  posYInitial=degradPos.getY();
    G4double  posZInitial=degradPos.getZ();
    G4double  timeInitial=degradTime;
    G4String line;
    std::vector<G4double> v;
    
    std::ifstream inFile;
    G4String fname= "DEGRAD.OUT";
    inFile.open(fname,std::ifstream::in);
    
    G4cout<< "Working in "<<fname<<G4endl;
    
    nline=1;
    electronNumber=0;
    while (getline(inFile, line,'\n')) {
        
        std::istringstream iss(line);//stream de strings
        
        if (nline ==1) {
            while (iss >> n) {
                v.push_back(n);
            }
            
            eventNumber=v[0];
            Nep=v[1];
            //  Nexc=v[2];
            v.clear();
        }
        //Ionizations
        if (nline ==2)  {
            
            fastStep.SetNumberOfSecondaryTracks(1E6); // reasonable max # of electrons created by degrad
            
            while (iss >> n) {
                v.push_back(n); //o n é adicionado ao vector
            }
            
            for (i=0;i<v.size();i=i+7){
                posXDegrad=v[i];
                posYDegrad=v[i+1];
                posZDegrad=v[i+2];
                timeDegrad=v[i+3];
                //convert from um to mm in GEANT4
                //also Y and Z axes are swapped in GEANT4 and Garfield++ relatively to Degrad
                posX=posXDegrad*0.001+posXInitial;
                posY=posZDegrad*0.001+posYInitial;
                posZ=posYDegrad*0.001+posZInitial;
                //std::cout << "DegradModel::DoIt(): v[i-4]" << v[i] << "," << v[i+1] << "," << v[i+2] << "," << v[i+3] << "," << v[i+4]   << std::endl;
                //std::cout << "DegradModel::DoIt(): xinitial, poxXDegrad [mm]" << posXInitial << ", " << posXDegrad*0.001 << std::endl;
    
                //convert ps to ns
                time=timeDegrad*0.001+timeInitial;
                
                
                G4ThreeVector myPoint;
                myPoint.setX(posX);
                myPoint.setY(posY);
                myPoint.setZ(posZ);
                
                //Check in which Physical volume the point bellongs
                G4Navigator* theNavigator= G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
                
                G4VPhysicalVolume* myVolume = theNavigator->LocateGlobalPointAndSetup(myPoint);
              
                G4String solidName=myVolume->GetName();
                
                if (G4StrUtil::contains(solidName,"FIELDCAGE") || G4StrUtil::contains(solidName,"GAS") ){

                    electronNumber++;
                    XenonHit* xh = new XenonHit();
                    xh->SetPos(myPoint);
                    xh->SetTime(time);
                    fGasBoxSD->InsertXenonHit(xh);
                    
                    // Create secondary electron
                    //if(electronNumber % 50 == 0){     // comment this condition, EC, 2-Dec-2021.
                    // G4DynamicParticle electron(G4Electron::ElectronDefinition(),G4RandomDirection(), 7.0*eV);

                    G4DynamicParticle electron(NEST::NESTThermalElectron::ThermalElectronDefinition(),G4RandomDirection(), 1.13*eV);
                    G4Track *newTrack=fastStep.CreateSecondaryTrack(electron, myPoint, time,false);

                    // }
                }
            }
            v.clear(); //Faz reset ao vector senão vai continuar a adicionar os dadosadicionar os dados
            nline=0;
            
        }
        nline++;
        
        
    }
    inFile.close();
    G4cout << "Number of initial electrons: " << electronNumber << G4endl;
    
    
}


