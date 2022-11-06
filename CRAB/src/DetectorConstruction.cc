#include "DetectorConstruction.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
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
#include "DegradModel.hh"
#include "GarfieldVUVPhotonModel.hh"
#include "G4SDManager.hh"


DetectorConstruction::DetectorConstruction(GasModelParameters* gmp)
    :
    fGasModelParameters(gmp),
    checkOverlaps(1),
    worldHalfLength(3.*m), //World volume is a cube with side length = 3m;
    wallThickness(0.05*m), //thickness of the aluminum walls
    caloThickness(1.*mm), // thickness of the silicon detectors
    gasPressure(1.*bar), // Pressure inside the gas
    temperature(273.15*kelvin), // temperature
    neonPercentage(85.72), // mixture settings
    co2Percentage(9.52)
{
  detectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction() {
  delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
  /* The World volume is a vacuum in which a gastube is placed with the walls made out of Aluminum. The
  endcaps are Silicon detectors, used as calorimeter 
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
  const G4double torr = 1. / 760. * atmosphere;
  std::cout << "DetectorConstruction::Construct() pressure [bar,pascals?] is " << gasPressure/bar << "," << gasPressure << std::endl;
  //  detectorMessenger->

  //Gas material: mixture of HeIso or ArCO2
  G4double nMoles = gasPressure/torr / (8.314 * joule / mole * temperature);
  G4Material* mixture=NULL;
  G4VPhysicalVolume* physiWorld = NULL;
  
  G4Material* air = man->FindOrBuildMaterial("G4_AIR");
  G4Material* lead = man->FindOrBuildMaterial("G4_Pb");
  //  const static G4double Torr = gasPressure/torr* 1. / 760. * atmosphere;
  G4Material* Xenon = man->ConstructNewGasMaterial ("XenonGas", "G4_Xe", 296.*kelvin, gasPressure, false);
  G4Material* glass= man->FindOrBuildMaterial("G4_MAGNESIUM_FLUORIDE");
  G4Material* kapton=man->FindOrBuildMaterial("G4_KAPTON");
  G4Material* fSteel = new G4Material("StainlessSteel", 7.80 * g/cm3, 3);
  G4Element* elFe = man->FindOrBuildElement("Fe");
  G4Element* elNi = man->FindOrBuildElement("Ni");
  G4Element* elCr = man->FindOrBuildElement("Cr");
  fSteel->AddElement(elFe, 70 * perCent);
  fSteel->AddElement(elCr, 18 * perCent);
  fSteel->AddElement(elNi, 12 * perCent);
  G4Element* elAl = man->FindOrBuildElement("Al");
  G4Element* elO = man->FindOrBuildElement("O");
  G4Material* fMacor=new G4Material("ceramic", 2.52 * g/cm3, 2);
  fMacor->AddElement(elAl,2);
  fMacor->AddElement(elO,3);
  
  G4Material* steel304=man->FindOrBuildMaterial("StainlessSteel");
  G4Material* ceramic=man->FindOrBuildMaterial("ceramic");  //MACOR

  const G4int nEntriesXenonIndex = 11;
  G4double photonEnergyXenonIndex[nEntriesXenonIndex]={6.25*eV,6.41*eV,6.58*eV,6.75*eV,6.93*eV,7.12*eV,7.32*eV,7.54*eV,7.77*eV,8.01*eV,8.27*eV};
  
  G4double XenonRindex[nEntriesXenonIndex]={1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
  
  G4MaterialPropertiesTable* mptXenon= new G4MaterialPropertiesTable();
  mptXenon->AddProperty("RINDEX",photonEnergyXenonIndex,XenonRindex,nEntriesXenonIndex);
  Xenon->SetMaterialPropertiesTable(mptXenon);
  
  //CAM GLASS Refractive index
  const G4int nEntriesMgF2Index = 9;
  G4double photonEnergyMgF2[nEntriesMgF2Index]={6.19*eV,6.44*eV,6.70*eV,6.97*eV,7.26*eV,7.55*eV,7.86*eV,8.18*eV,8.51*eV};
  

  G4double MgF2Rindex[nEntriesMgF2Index]={1.423,1.428,1.434,1.440,1.447,1.455,1.465,1.476,1.489};
  G4double MgF2Reflectivity[nEntriesMgF2Index]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  G4double MgF2Efficiency[nEntriesMgF2Index]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  
  
  G4MaterialPropertiesTable* mptGlass= new G4MaterialPropertiesTable();
  mptGlass->AddProperty("RINDEX",photonEnergyMgF2,MgF2Rindex,nEntriesMgF2Index);
  mptGlass-> AddProperty("REFLECTIVITY",photonEnergyMgF2,MgF2Reflectivity,nEntriesMgF2Index);
  mptGlass-> AddProperty("EFFICIENCY",photonEnergyMgF2,MgF2Efficiency,nEntriesMgF2Index);

  glass->SetMaterialPropertiesTable(mptGlass);
  
  
  
  //--------------------------------------------------
  // CAM PHOTOCATHODE
  //--------------------------------------------------
  
  
  G4double photonEnergyCAM[]={6.26*eV,6.31*eV,6.38*eV,6.42*eV,6.45*eV,6.49*eV,6.52*eV,6.55*eV,6.59*eV,6.61*eV,6.66*eV,6.70*eV,6.75*eV,6.81*eV,6.85*eV,6.89*eV,6.95*eV,6.98*eV,7.02*eV,7.08*eV,7.15*eV,7.21*eV,7.30*eV,7.37*eV,7.43*eV,7.50*eV,7.56*eV,7.60*eV,7.68*eV,7.73*eV,7.77*eV,7.86*eV,7.95*eV,8.04*eV,8.10*eV,8.14*eV};
  
  
  
  
  G4double photocath_EFF[]=
  {3.96,5.20,6.63,7.96,8.99,9.85,10.46,11.46,12.18,12.94,14.17,14.61,16.00,17.01,18.07,18.63,19.79,20.40,21.68,21.68,23.04,23.04,23.75,25.24,26.01,26.01,27.64,28.49,29.37,29.37,30.28,32.17,33.17,35.24,36.33,37.45}; //Enables 'detection' of photons
  
  assert(sizeof(photocath_EFF) == sizeof(photonEnergyCAM));
  
  
  G4double reflectivityPhotocathode[] =
  { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  assert(sizeof(reflectivityPhotocathode) == sizeof(photonEnergyCAM));
  
  

  gasboxR = 15.*cm;
  gasboxH = 60.*cm;

  G4double halfColimatorLength=0.2*cm;
  G4double rIntColimator=1.0*mm;
  G4double rIntColimator2=3.75*mm;
  G4double rExtColimator=5.*cm;
  //detectorHalfZ=5*cm;
  G4double camHalfLength=0.25*cm;
  G4double camRadius=2.4*cm;
  G4double distColimattorToDetector=2.5*cm;
  G4double macorHalfY=1.3*cm/2;
  G4double windowHalfLength=25.*um;
  
  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  G4VSolid* worldSolid
  = new G4Box("worldBox",500.*cm,400.*cm,500.*cm);
  G4LogicalVolume* worldLogical  = new G4LogicalVolume(worldSolid,air,"worldLogical");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                      false,0,checkOverlaps);
  
  
  // collimator Tube
  G4VSolid* collimatorSolid = new G4Tubs("collimatorTube",rIntColimator,rExtColimator,halfColimatorLength,0.,twopi);
  G4LogicalVolume* collimatorLogical = new G4LogicalVolume(collimatorSolid,lead,"collimatorLogical");
  
  
  
  G4RotationMatrix rotm  = G4RotationMatrix();
  rotm.rotateX(90*deg);
  rotm.rotateY(0*deg);
  rotm.rotateZ(0*deg);
  G4RotationMatrix rotnull  = G4RotationMatrix();
  rotm.rotateX(0*deg);
  rotm.rotateY(0*deg);
  rotm.rotateZ(0*deg);
  
  G4ThreeVector position = G4ThreeVector(0.,-halfColimatorLength-distColimattorToDetector,0.*mm);
  
  G4Transform3D transform = G4Transform3D(rotm,position);
  
  
  /* 
  new G4PVPlacement(transform,collimatorLogical,"collimatorPhysical",worldLogical,
                                                            false,0,checkOverlaps);
  */
  
  
  //Make XenonVolume
  
  G4VSolid* detectorSolid=new  G4Tubs("detectorBox",0,gasboxR,gasboxH*0.5,0.,twopi);
  
  ////Place Detector in world
  logicGasBox = new G4LogicalVolume(detectorSolid,Xenon,"detectorLogical");
  
  G4ThreeVector position3 = G4ThreeVector(0.,gasboxH*0.5,0.*mm);
  G4Transform3D transform3 = G4Transform3D(rotm,position3);
  
  G4VPhysicalVolume* detectorPhysical=    new G4PVPlacement(transform3,             //rotation and position
                                                            logicGasBox,            //its logical volume
                                                            "detectorPhysical",             //its name
                                                            worldLogical,             //its mother  volume
                                                            true,                 //no boolean operation
                                                            0,                 //copy number
                                                            checkOverlaps);       // checking overlaps
  
  
  //Macor
  
  
  G4VSolid* macorSolid=new  G4Tubs("macorTube",gasboxR,gasboxR+rExtColimator,macorHalfY,0.,twopi);
  
  ////Place macor in world
  G4LogicalVolume* macorLogical = new G4LogicalVolume(macorSolid,ceramic,"macorLogical");
  
  G4ThreeVector position4 = G4ThreeVector(0.,macorHalfY,0.*mm);
  G4Transform3D transform4 = G4Transform3D(rotm,position4);
  /*
  new G4PVPlacement(transform4,             //rotation and position
                    macorLogical,            //its logical volume
                    "macorPhysical",             //its name
                    worldLogical,             //its mother  volume
                    true,                 //no boolean operation
                    0,                 //copy number
                    checkOverlaps);       // checking overlaps
  */
  //CAM
  G4VSolid* camSolid = new G4Tubs("camWindow",0.,camRadius,camHalfLength,0.,twopi);
  G4LogicalVolume* camLogical = new G4LogicalVolume(camSolid,glass,"camLogical");
  
  
  //Make camLogical mother and photocathode daughter
  G4ThreeVector positionCAM = G4ThreeVector(gasboxR/2.0,0.*mm,+gasboxH*0.5-0.5*camHalfLength); // just inside the volume
  G4Transform3D transformCAM = G4Transform3D(rotnull,positionCAM);
  G4VPhysicalVolume* camPhysical= new G4PVPlacement(transformCAM,camLogical,"camPhysical",logicGasBox,
                                                    false,0,checkOverlaps);  
  
  
  // lensTube
  const G4double lensTubeLength(10.0*cm);
  const G4double tubeThickness(1.0*mm);
  G4VSolid* lensTube = new G4Tubs("lensTube",camRadius-tubeThickness,camRadius,lensTubeLength/2.,0.,twopi);
  G4LogicalVolume* lensTubeLogical = new G4LogicalVolume(lensTube,steel304,"lensTube");
  G4ThreeVector positionLT = G4ThreeVector(gasboxR/2.0,0.0,+gasboxH*0.5-lensTubeLength/2.-camHalfLength*2);  
  G4Transform3D transformLT = G4Transform3D(rotnull,positionLT);
  //  G4VPhysicalVolume* lensTubePhysical = new G4PVPlacement(transformLT,lensTubeLogical,"lensTubeLogical",logicGasBox,
  //                                                               false,0,checkOverlaps);
  // lens
  const G4double lensRcurve (4*camRadius); // 
  //  const G4double lensThickness(6.9*mm);
  //  const G4double lensRadius (10.0*mm);
  G4double s2(lensRcurve); // approximately
  G4double nind(MgF2Rindex[2]); // pick an index of refrac that's typical
  G4double s1;
  s1 = gasboxH -0.5*camHalfLength - s2;
  G4double f(1.2 * lensRcurve/(2.0*(nind-1)));  // 1.2 is empirical in order to force focussing of parallel rays. EC, 16-Aug-2022.
  s2 = 1.0/(1/f-1.0/s1);
  std::cout<<"s1,s2 [mm]: "<< s1 << ", " << s2 << std::endl;
  
  const G4ThreeVector posLensTubeIntersect (0.,lensRcurve*1.5,0.);
  G4Orb* sLensTube = new G4Orb("sLensSphereTube",lensRcurve); //new G4Tubs("sLensTube",0.0,lensRadius,lensThickness*2.0,0.0,twopi);
  G4Orb* sLensOrb = new G4Orb("sLensSphere",lensRcurve);
  G4IntersectionSolid* sLens =  new G4IntersectionSolid("sLens",sLensTube,sLensOrb,&rotnull,posLensTubeIntersect);
  G4LogicalVolume* lensLogical = new G4LogicalVolume(sLens,glass,"Lens");
  G4ThreeVector positionLens = G4ThreeVector(gasboxR/2.0,0.,+gasboxH*0.5-s2-camHalfLength*2);  
  G4Transform3D transformLens = G4Transform3D(rotm,positionLens);
  G4VPhysicalVolume* lensPhysical = new G4PVPlacement(transformLens,lensLogical,"Lens",logicGasBox,false,0,checkOverlaps);


  //OPTICAL SURFACES
  
  //Xenon-GLASS

    // commenting out does not suppress reflections. EC, 12-Apr-2022.
  
  G4OpticalSurface* opXenon_Glass = new G4OpticalSurface("XenonCamSurface");
  opXenon_Glass->SetModel(glisur);                  // SetModel
  opXenon_Glass->SetType(dielectric_dielectric);   // SetType
  opXenon_Glass->SetFinish(ground);                 // SetFinish
  opXenon_Glass->SetPolish(0.0);
  new G4LogicalBorderSurface("XenonCamSurface",detectorPhysical,camPhysical,opXenon_Glass);

  G4OpticalSurface* opXenon_Glass2 = new G4OpticalSurface("XenonLensSurface");
  opXenon_Glass2->SetModel(glisur);                  // SetModel
  opXenon_Glass2->SetType(dielectric_dielectric);   // SetType
  opXenon_Glass2->SetFinish(polished);                 // SetFinish
  opXenon_Glass2->SetPolish(0.0);
  new G4LogicalBorderSurface("XenonLensSurface",detectorPhysical,lensPhysical,opXenon_Glass2);
  //    new G4LogicalSkinSurface("XenonLensSurface",lensLogical,opXenon_Glass2);

  
  //    // visualization attributes ------------------------------------------------
  
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  //  collimatorLogical->SetVisAttributes(red);
  //  collimatorLogical2->SetVisAttributes(green);
  logicGasBox->SetVisAttributes(blue);
  camLogical->SetVisAttributes(yellow);
  lensTubeLogical->SetVisAttributes(purple);
  lensLogical->SetVisAttributes(purple);


  //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
  G4Region* regionGas = new G4Region("GasRegion");
  regionGas->AddRootLogicalVolume(logicGasBox);
    
  return physiWorld;

}

void DetectorConstruction::ConstructSDandField(){
  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  G4String GasBoxSDname = "interface/GasBoxSD";
  GasBoxSD* myGasBoxSD = new GasBoxSD(GasBoxSDname);
  SDManager->SetVerboseLevel(1);
  SDManager->AddNewDetector(myGasBoxSD);
  SetSensitiveDetector(logicGasBox,myGasBoxSD);
  
  //These commands generate the four gas models and connect it to the GasRegion
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
  new DegradModel(fGasModelParameters,"DegradModel",region,this,myGasBoxSD);
  new GarfieldVUVPhotonModel(fGasModelParameters,"GarfieldVUVPhotonModel",region,this,myGasBoxSD);

}

