#include "DetectorConstruction.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4OpticalSurface.hh"
#include "G4Trd.hh"
#include "G4Threading.hh"
#include "G4RegionStore.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Cons.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trd.hh"
#include "DetectorMessenger.hh"
#include "CRAB_HybridGeo.h"
#include "GasBoxSD.hh"
#include "DegradModel.hh"
#include "GarfieldVUVPhotonModel.hh"
#include "G4SDManager.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4MultiUnion.hh"
#include "Visibilities.hh"
#include "HexagonMeshTools.hh"

#include "SampleFromSurface.hh"
#include "config.h"

#ifdef With_Opticks
#include "G4CXOpticks.hh"
#include <cuda_runtime.h>
#include "SEventConfig.hh"
#endif

#include <chrono>

// Lets use tetgen
#ifdef WITH_SYS_TETGEN
#define USE_CADMESH_TETGEN 1
#endif
#include "CADMesh.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidTube.hh"
#include "CRAB_CSG.hh"
#include "SimpleBoxGeo.hh"
#include "SensorSD.hh"

DetectorConstruction::DetectorConstruction(GasModelParameters* gmp) :
    fGasModelParameters(gmp),
    Lab_size(1 *m),
    chamber_diam(16.4 * cm),
    chamber_length (43.18* cm), // Config files vary
    chamber_thickn (7. * mm),
    gas_pressure_(10 * bar),
    FielCageGap(21.26*cm),
    Active_diam(8.6 * cm),
    temperature(293*kelvin)
{
    detectorMessenger = new DetectorMessenger(this);

}

DetectorConstruction::~DetectorConstruction() {
    delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    //Materials

    G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());

    // Constructing Lab Space
    G4String lab_name="LAB";
    G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2,Lab_size/2,Lab_size/2);



    G4LogicalVolume * lab_logic_volume = new G4LogicalVolume(lab_solid_volume,air,lab_name) ; //G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR")

    G4double Offset=0.2*cm;
    // Placement of the Items
    auto labPhysical = new G4PVPlacement(0, G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0, false,0,false);

    // Visuals


#ifdef With_Opticks
    CRAB_CSG *geo=new CRAB_CSG(fGasModelParameters);
    //SimpleBoxGeo *geo = new SimpleBoxGeo();
#endif

#ifndef With_Opticks
    CRAB_HybridGeo *geo=new CRAB_HybridGeo(fGasModelParameters);
#endif

    geo->SetMotherLab(lab_logic_volume);
    geo->SetOffset(0.2*cm);
    //geo->SetYield(10);
    geo->SetGasPressure(gas_pressure_);
    geo->SetTemperature(temperature);
    geo->Construct();
    gas_logic=geo->GasLogic();
    //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
    G4Region* regionGas = new G4Region("GasRegion");
    regionGas->AddRootLogicalVolume(gas_logic);





#ifdef With_Opticks
    std::cout <<"Setting our detector geometry with opticks" <<std::endl;
    G4CXOpticks::Get()->SetGeometry(labPhysical);
    std::cout << SEventConfig::Desc() <<std::endl;
#endif

    return labPhysical;

}
void DetectorConstruction::ConstructSDandField(){

    G4SDManager* SDManager = G4SDManager::GetSDMpointer();
    G4String GasBoxSDname = "interface/GasBoxSD";

    GasBoxSD* myGasBoxSD = new GasBoxSD(GasBoxSDname);
    SDManager->SetVerboseLevel(1);
    SDManager->AddNewDetector(myGasBoxSD);
    SetSensitiveDetector(gas_logic,myGasBoxSD);


    //These commands generate the four gas models and connect it to the GasRegion
    G4Region* region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
    new DegradModel(fGasModelParameters,"DegradModel",region,this,myGasBoxSD);
    new GarfieldVUVPhotonModel(fGasModelParameters,"GarfieldVUVPhotonModel",region,this,myGasBoxSD);
}

