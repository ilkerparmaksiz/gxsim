//
// Created by Ilker Parmaksiz on 2/15/24.
//

#include "SimpleBoxGeo.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "config.h"
#ifdef With_Opticks
#include <cuda_runtime.h>
#include "SEventConfig.hh"
#include "G4CXOpticks.hh"
#include "U4Physics.hh"
#include <cuda_runtime.h>
#include <globals.hh>
#include "U4Scint.h"
#include "U4Material.hh"
#endif
#include "MaterialsList.hh"
#include "OpticalMaterialProperties.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "U4SensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
SimpleBoxGeo::SimpleBoxGeo() {

    // Initialise the values
    chamber_diam=16.4 * cm;
    chamber_length=43.18* cm;
    gas_pressure_=10 * bar;
    FielCageGap=21.26*cm;
    Active_diam=8.6 * cm;
}
SimpleBoxGeo::~SimpleBoxGeo(){}
void SimpleBoxGeo::Construct() {
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4Material *GXenon = materials::GXe(1 * bar, 293 * kelvin);
    GXenon->SetMaterialPropertiesTable(opticalprops::GXe(GetGasPressure(), GetTemperature(), GetYield()));

    G4Material *MgF2 = materials::MgF2();
    MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());


    G4Material *Steel = materials::Steel();
    Steel->SetMaterialPropertiesTable(opticalprops::STEEL());

    G4double env_sizeXY = 20 * cm, env_sizeZ = 30 * cm;
    G4double world_sizeXY = 1.2 * env_sizeXY;
    G4double world_sizeZ = 1.2 * env_sizeZ;



    auto CubeDetector_Solid = new G4Box("DetectorSolid", 1 * cm, 1 * cm, 1* mm);  // its size
    auto CubeDetectorLogic = new G4LogicalVolume(CubeDetector_Solid, MgF2, "Detector_Logic");



    G4LogicalVolume * logicWorld=Mother;

    //
    // Envelope
    //
    auto solidEnv = new G4Box("GasXe_Solid",                    // its name
                              0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

    auto logicEnv = new G4LogicalVolume(solidEnv, GXenon, "GAS_");


    auto Steel_Cover_solid = new G4Box("SteelCover_Solid",                    // its name
                                       0.5 * env_sizeXY + 5 * CLHEP::mm, 0.5 * env_sizeXY + 10 * CLHEP::mm,
                                       0.5 * env_sizeZ + 5 * CLHEP::mm);  // its size
    auto Steel_Cover_logic = new G4LogicalVolume(Steel_Cover_solid,  // its solid
                                                 Steel,                                     // its material
                                                 "SteelCover_logic");                                 // its name
    //
    //always return the physical World
    //
    auto SteelPlace = new G4PVPlacement(0, G4ThreeVector(0, 0, 0) ,Steel_Cover_logic,"SteelCover", Mother, 0, 0,0);
    auto GasXePlace = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GasXe", logicEnv, SteelPlace, 0, 0, 0);
    auto DetectorPlace = new G4PVPlacement(0, G4ThreeVector(0, 0, -10 * cm), "Detector", CubeDetectorLogic, GasXePlace, 0, 0,0);

    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
    Steel_Cover_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VisAttributes gasVis = new G4VisAttributes();
    gasVis.SetColor(0, 1, 1);
    gasVis.SetForceCloud(true);
    logicEnv->SetVisAttributes(gasVis);
    G4VisAttributes detVis = new G4VisAttributes();
    detVis.SetColor(0, 0, 1);
    detVis.SetForceCloud(true);
    CubeDetectorLogic->SetVisAttributes(detVis);
    //G4OpticalSurface *OpSteelSurf = new G4OpticalSurface("SteelSurface", unified, polished, dielectric_metal);
    //OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());


    G4OpticalSurface *opXenon_Glass = new G4OpticalSurface("DetectorSurface");
    opXenon_Glass->SetMaterialPropertiesTable(opticalprops::MgF2());
    opXenon_Glass->SetModel(glisur);                  // SetModel
    opXenon_Glass->SetType(dielectric_dielectric);   // SetType
    opXenon_Glass->SetFinish(ground);                 // SetFinish
    opXenon_Glass->SetPolish(0);

    new G4LogicalBorderSurface("DetectorSurface", GasXePlace, DetectorPlace, opXenon_Glass);

    Detector=CubeDetectorLogic;
    G4SDManager *Mang = G4SDManager::GetSDMpointer();
    //PhotonSD *SD = new PhotonSD("PMT1");
    //Mang->AddNewDetector(SD);
    //Detector->SetSensitiveDetector(SD);
    this->SetGasLogic(logicEnv);


}