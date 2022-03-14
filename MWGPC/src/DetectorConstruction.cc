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
    checkOverlaps(1),
    worldHalfLength(3.*m), //World volume is a cube with side length = 3m;
    wallThickness(0.002*m), //thickness of the Copper walls
    gasPressure(10.*bar), // Pressure inside the gas ..... also used for magboltz, EC 28-Oct-2021.
    numHexes(1), // number of anode wires
    temperature(273.15*kelvin), // temperature
    argonPercentage(90.00), // mixture settings
    ch4Percentage(10.00)
{
  detectorMessenger = new DetectorMessenger(this);
  fHDEt = 0;

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

  gasPressure = GetGasPressure();
  temperature = GetTemperature();
  
  //Gas material: P10
  G4double nMoles = gasPressure / (8.314 * joule / mole * temperature);
  G4Material* mixture=NULL;
  G4VPhysicalVolume* physiWorld = NULL;
  // TPC setup
  gasboxR = 0.0325*m;
  gasboxH = 0.19*m;
  G4double wireR = 80 * um;
  
  G4Element* elC = man->FindOrBuildElement("C");
  G4Element* elH = man->FindOrBuildElement("H");
  G4Element* elAr = man->FindOrBuildElement("Ar");
  
  G4double molarMass = 40.*g/mole;  // pure Ar
  
  G4double gasDensityAr = nMoles * molarMass;
  G4cout << "gasPressure: " << G4BestUnit(gasPressure, "Pressure")
     << G4endl;
  G4cout << "gasTemperature: " << G4BestUnit(temperature, "Temperature")
     << G4endl;
  G4cout << "numHexes: " << std::to_string(numHexes) 
     << G4endl;

  G4Material* Argon = new G4Material("ar", 10, molarMass, gasDensityAr,
                                      kStateGas, temperature, gasPressure);
  G4double molfracAr = (argonPercentage/100.) * molarMass;

  // CH4
  molarMass = 16.0*g/mole;  // source wikipedia
  G4double gasDensityCH4 = nMoles * molarMass;
  G4Material* CH4 = new G4Material("ch4", gasDensityCH4, 2,
                                    kStateGas, temperature, gasPressure);
  CH4->AddElement(elC, 1);
  CH4->AddElement(elH, 4);
  
  G4double molfracCH4 = (ch4Percentage/100.)*molarMass;
  

  G4double molfracAr_norm = molfracAr/(molfracAr+molfracCH4);
  G4double molfracCH4_norm = molfracCH4/(molfracAr+molfracCH4);


  G4double gasDensityMixture = (argonPercentage/100.) * gasDensityAr +
    ch4Percentage/100. * gasDensityCH4;



  // Coment in one or the other code chunk below. EC, 26-Oct-2021.
  // Also remember to have correct gasfile name in *.mac. And code commented in HeedModel.cc if not creating the gasfile.

  
  /*
  G4cout << "gasDensityAr: " << G4BestUnit(gasDensityAr, "Volumic Mass")
     << G4endl;
  G4cout << "gasDensityCH4: " << G4BestUnit(gasDensityCH4, "Volumic Mass")
         << G4endl;
  G4cout << "Molar fraction Ar: " << molfracAr_norm << G4endl;
  G4cout << "Molar fraction CH4: " << molfracCH4_norm << G4endl;


  mixture = new G4Material("mixture", gasDensityMixture, 2);
  
  mixture->AddMaterial(Argon, molfracAr_norm);
  mixture->AddMaterial(CH4, molfracCH4_norm);

  G4cout << "gasDensityP10: " << G4BestUnit(gasDensityMixture, "Volumic Mass") << G4endl;
  */
  
  
  
  G4String symbol, name;
  G4double density = 1.17*mg/cm3 * gasPressure/bar;
  std::cout << "DetectorConstruction: density is " << G4BestUnit(density, "Volumic Mass")  << std::endl;
  G4int ncomponents;

  //  mixture =
  //  new G4Material(name="mixture",density,ncomponents=2);//,kStateGas, temperature,gasPressure);
  mixture = new G4Material("mixture", density, 2);
  G4Material* O2 =  man->FindOrBuildMaterial("G4_O");
  G4Material* N2 =  man->FindOrBuildMaterial("G4_N");
  mixture->AddMaterial(N2, 0.99);
  mixture->AddMaterial(O2, 0.01);
  //    mixture->AddMaterial(N2, 1.00);
  
  
  
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
  std::cout << "DetCon:: gas radius, half-height [mm]: " <<gasboxR << ", " << gasboxH*0.5 << std::endl;
  G4Tubs* solidGasBox = new G4Tubs("solid_gasbox_tube",0.0/*wireR*/,gasboxR,gasboxH*0.5, 0., twopi);
  logicGasBox =
    new G4LogicalVolume(solidGasBox, mixture, "solidGasBox_log");
  new G4PVPlacement(myRotation,G4ThreeVector(), logicGasBox,"solidGasBox_phys",logicWorld,false,0,checkOverlaps);
  
  
  //Copper Wall
  G4Tubs* solidWalls = new G4Tubs("solid_tube_wall",gasboxR,gasboxR+wallThickness, gasboxH*0.5+wallThickness, 0., twopi);
  G4LogicalVolume* logicWall =
    new G4LogicalVolume(solidWalls, copperMaterial, "solidWall_log");
  new G4PVPlacement(myRotation,G4ThreeVector(), logicWall,
                    "solidWall_phys",logicWorld,false,0,checkOverlaps);
  logicGasBox->SetVisAttributes(blue);
  logicWall->SetVisAttributes(green);

  //Copper S End Wall
  G4Tubs* solidWallsS = new G4Tubs("solid_tube_wallS",0.,gasboxR+wallThickness,0.5*wallThickness, 0., twopi);
  G4LogicalVolume* logicWallS =
    new G4LogicalVolume(solidWallsS, copperMaterial, "solidWallS_log");
  new G4PVPlacement(myRotation,G4ThreeVector(0.,gasboxH*0.5+wallThickness*0.5,0.), logicWallS,
                    "solidWallS_phys",logicWorld,false,0,checkOverlaps);
  //logicGasBox->SetVisAttributes(blue);
  logicWallS->SetVisAttributes(red);

  //Copper N End Wall
  G4Tubs* solidWallsN = new G4Tubs("solid_tube_wallN",0.,gasboxR+wallThickness,0.5*wallThickness, 0., twopi);
  G4LogicalVolume* logicWallN =
    new G4LogicalVolume(solidWallsN, copperMaterial, "solidWallN_log");
  new G4PVPlacement(myRotation,G4ThreeVector(0.,-gasboxH*0.5-0.5-wallThickness*0.5,0.), logicWallN,
                    "solidWallN_phys",logicWorld,false,0,checkOverlaps);
  //logicGasBox->SetVisAttributes(blue);
  logicWallN->SetVisAttributes(red);

  
  //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimrghasulationModel
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


  fHDEt = new HeedDeltaElectronModel(fGasModelParameters,"HeedDeltaElectronModel",region,this,myGasBoxSD);
  // at this point fbins and fNumbins have values in the HeedModel track instance. Let us hold them in this detcon object rather than
  // needing to get them at EndEventAction, by which time trackstack is empty and these are unavailable.

  fBinSz = fHDEt->GetBinsz();
  fNumBins = fHDEt->GetNumbins();

  std::cout << "DetCon::ConSDandField(): this, binsz, numbins" << this << ", " << fBinSz << "," << fNumBins << std::endl;
}

