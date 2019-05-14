#include "DetectorConstruction.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalBorderSurface.hh"
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
#include "GasBoxSD.hh"
#include "SiliconSD.hh"
#include "HeedDeltaElectronModel.hh"
#include "HeedNewTrackModel.hh"
#include "G4SDManager.hh"


DetectorConstruction::DetectorConstruction(GasModelParameters* gmp)
    :
    fGasModelParameters(gmp),
    checkOverlaps(0),
    worldHalfLength(3.*m), //World volume is a cube with side length = 3m;
    wallThickness(0.002*m), //thickness of the Copper walls
    gasPressure(10.*bar), // Pressure inside the gas
    numHexes(1), // number of anode wires
    temperature(273.15*kelvin), // temperature
    argonPercentage(90.00), // mixture settings
    ch4Percentage(10.00)
{
  detectorMessenger = new DetectorMessenger(this);


}

DetectorConstruction::~DetectorConstruction() {
  delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
  /* The World volume is a vacuum in which a gastube is placed with the walls made out of copper. 
  */
    
  //Colors for visualization
  G4VisAttributes* red = new G4VisAttributes(G4Colour(1., 0., 0.));
  G4VisAttributes* green = new G4VisAttributes(G4Colour(0., 1., 0.));
  G4VisAttributes* blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  G4VisAttributes* yellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.));
  G4VisAttributes* purple = new G4VisAttributes(G4Colour(1.0, 0., 1.0));

  /*First: build materials
    World: vacuum
    Walls: Aluminum
    Gas: mixture of Helium and Isobutane or Ar and CO2
    Calorimeter: Silicon 
  */
  
    
  //World material: vacuum
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  G4Material* vacuum = man->FindOrBuildMaterial("G4_Galactic");
  
  //Gas material: P10
  G4double nMoles = gasPressure / (8.314 * joule / mole * temperature);
  G4Material* mixture=NULL;
  G4VPhysicalVolume* physiWorld = NULL;
  // TPC setup
  gasboxR = 0.0325*m;
  gasboxH = 0.19*m;

  G4Element* elC = man->FindOrBuildElement("C");
  G4Element* elH = man->FindOrBuildElement("H");
  G4Element* elAr = man->FindOrBuildElement("Ar");
  
  G4double molarMass = 40.*g/mole;  // pure Ar
  
  G4double gasDensityAr = nMoles * molarMass;
  G4cout << "gasPressure: " << G4BestUnit(gasPressure, "Pressure")
     << G4endl;
  G4cout << "gasDensityAr: " << G4BestUnit(gasDensityAr, "Volumic Mass")
     << G4endl;
  G4cout << "numHexes: " << std::to_string(numHexes) 
     << G4endl;

  G4Material* Argon = new G4Material("ar", 10, molarMass, gasDensityAr,
                                      kStateGas, temperature, gasPressure);
  G4double molfracAr = (argonPercentage/100.) * molarMass;

  // CH4
  molarMass = 16.0*g/mole;  // source wikipedia
  G4double gasDensityCH4 = nMoles * molarMass;
  G4cout << "gasDensityCH4: " << G4BestUnit(gasDensityCH4, "Volumic Mass")
         << G4endl;
  G4Material* CH4 = new G4Material("ch4", gasDensityCH4, 2,
                                    kStateGas, temperature, gasPressure);
  CH4->AddElement(elC, 1);
  CH4->AddElement(elH, 4);
  
  G4double molfracCH4 = (ch4Percentage/100.)*molarMass;
  

  G4double molfracAr_norm = molfracAr/(molfracAr+molfracCH4);
  G4double molfracCH4_norm = molfracCH4/(molfracAr+molfracCH4);

  G4cout << "Molar fraction Ar: " << molfracAr_norm << G4endl;
  G4cout << "Molar fraction CH4: " << molfracCH4_norm << G4endl;


  G4double gasDensityMixture = (argonPercentage/100.) * gasDensityAr +
    ch4Percentage/100. * gasDensityCH4;
  
  mixture = new G4Material("mixture", gasDensityMixture, 2);
  
  mixture->AddMaterial(Argon, molfracAr_norm);
  mixture->AddMaterial(CH4, molfracCH4_norm);

  G4cout << "gasDensityP10: " << G4BestUnit(gasDensityMixture,
                                                "Volumic Mass") << G4endl;
  
  //geometry dimensions:
  //Aluminum walls
  G4Material* copperMaterial = man->FindOrBuildMaterial("G4_Cu");
  
  //World Volume
  G4Box* solidWorld = new G4Box("solidWorld_box", worldHalfLength, worldHalfLength, worldHalfLength);
  G4LogicalVolume* logicWorld =
  new G4LogicalVolume(solidWorld, vacuum, "solidWorld_log");
  
  physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,
                                                    "solidWorld_phys", 0, false, 0, checkOverlaps);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  

  // rotation makes cylinder axis lie along y axis.

  //GasBox volume
  G4RotationMatrix* myRotation = new G4RotationMatrix();
  myRotation->rotateX(90.*deg);
  myRotation->rotateY(0.*deg);
  myRotation->rotateZ(0.*rad);
  G4Tubs* solidGasBox = new G4Tubs("solid_gasbox_tube",0,gasboxR,gasboxH*0.5, 0., twopi);
  logicGasBox =
  new G4LogicalVolume(solidGasBox, mixture, "solidGasBox_log");
  new G4PVPlacement(myRotation,G4ThreeVector(), logicGasBox,"solidGasBox_phys",logicWorld,false,0,checkOverlaps);
  
  
  //Copper Wall
  G4Tubs* solidWalls = new G4Tubs("solid_tube_wall",gasboxR,gasboxR+wallThickness,gasboxH*0.5+wallThickness, 0., twopi);
  G4LogicalVolume* logicWall =
  new G4LogicalVolume(solidWalls, copperMaterial, "solidWall1_log");
  new G4PVPlacement(myRotation,G4ThreeVector(), logicWall,
                    "solidWall_phys",logicWorld,false,0,checkOverlaps);

  //logicGasBox->SetVisAttributes(blue);
  logicWall->SetVisAttributes(red);

  
  //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
  // THIS REGION IS THE ONE IN HeedModel.cc IN WHICH WE INSERT ENERGIZED WIRES AND MAGBOLTZ GAS.
  // NOTE THE call below in ConstructSDandField() to create local variable called region from this "GasRegion" string,
  // which is then passed to HeedNewTrackModel constructor.

  G4Region* regionGas = new G4Region("GasRegion");
  regionGas->AddRootLogicalVolume(logicGasBox);
    
  return physiWorld;

}


// This is a derived method from G4VUserDetectorconstruction, which is called by G4 somewhere.
void DetectorConstruction::ConstructSDandField(){
  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  G4String GasBoxSDname = "interface/GasBoxSD";
  GasBoxSD* myGasBoxSD = new GasBoxSD(GasBoxSDname);
  SDManager->SetVerboseLevel(1);
  SDManager->AddNewDetector(myGasBoxSD);
  SetSensitiveDetector(logicGasBox,myGasBoxSD);


  //These commands generate the four gas models and connect it to the GasRegion
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
  new HeedNewTrackModel(fGasModelParameters,"HeedNewTrackModel",region,this,myGasBoxSD);
  new HeedDeltaElectronModel(fGasModelParameters,"HeedDeltaElectronModel",region,this,myGasBoxSD);
}

