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
#include "OldCRAB.hh"


DetectorConstruction::DetectorConstruction(GasModelParameters* gmp) :
    fGasModelParameters(gmp),
    checkOverlaps(1),
    temperature(300*kelvin), // temperature
    Lab_size(5. *m),
    Gas_diam(19.4*cm),
    Gas_length(42.5*cm),
    chamber_diam   (16.4 * cm),
    chamber_length (43.18* cm), // Config files vary
    chamber_thickn (7. * mm),
    gas_pressure_(10 * bar),
    vtx_(0*cm,-1.6*cm,-5.25*cm),
    Active_diam(8.6 * cm),
    sc_yield_(25510./MeV),
    e_lifetime_(1000. * ms),
    HideSourceHolder_(false),
    ElGap_(7*mm),
    ELyield_(970/cm),
    fOffset(-0.8*cm)
{
    detectorMessenger = new DetectorMessenger(this);
    Sampler=std::make_shared<util::SampleFromSurface>(util::SampleFromSurface("Needles"));
}

DetectorConstruction::~DetectorConstruction() {
    delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    //Materials
    G4Material *gxe    = materials::GXe(gas_pressure_,68);
    G4Material *MgF2   = materials::MgF2();
    G4Material *Steel  = materials::Steel();
    G4Material *PEEK  = materials::PEEK();
    G4Material *vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material *teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    std::string crabpath= getenv("CRABPATH");

    // Constructing Lab Space
    G4String lab_name="LAB";
    G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2,Lab_size/2,Lab_size/2);

    //G4Tubs* gas_solid =new G4Tubs("GAS_", 0., chamber_diam/2., chamber_length/2. + chamber_thickn+3*cm/2, 0., twopi);
    G4Tubs* gas_solid =new G4Tubs("GAS_", 0., Gas_diam/2., Gas_length/2. + chamber_thickn, 0., twopi);


    // Optical Properties Assigned here
    MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());
    vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
    gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68,sc_yield_,e_lifetime_));

    G4LogicalVolume * lab_logic_volume = new G4LogicalVolume(lab_solid_volume,air,lab_name) ; //G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR")

    gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS_");

    // CAMERA WINDOW
    G4double camHalfLength=0.5*mm;
    G4double camRadius= 12.7*mm;
    G4VSolid* camSolid = new G4Tubs("camWindow",0.,camRadius,camHalfLength,0.,twopi);
    G4LogicalVolume* camLogical = new G4LogicalVolume(camSolid,MgF2,"Camera_logic");

    // lens
    const G4double lensRcurve (2.83*cm); // radius of curvature of MgF2 Lens
    const G4ThreeVector posLensTubeIntersect (0.,0.,-lensRcurve);
    // Create lens from the intersection of a sphere and a cylinder
    G4double maxLensLength = 4*mm;
    G4Tubs* sLensTube = new G4Tubs("sLensSphereTube", 0, 16.5 * mm/2, maxLensLength, 0.,twopi); // 4 mm is the max lens length
    G4Orb* sLensOrb = new G4Orb("sLensSphere",lensRcurve);
    G4IntersectionSolid* sLens =  new G4IntersectionSolid("sLens",sLensTube,sLensOrb, 0, posLensTubeIntersect);

    // Lens logical
    G4LogicalVolume* lensLogical = new G4LogicalVolume(sLens, MgF2, "Lens2");

    // Import Detector Geometry from an STL
    std::cout<< "CRAB Data " << crabpath <<std::endl;
   // auto FieldCage=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/FieldRings.stl");
    //auto Meshes=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Meshes.stl");
    auto Needle4=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_4cm.stl");
    auto Needle9=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_9cm.stl");
    auto Needle14=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_14cm.stl");
    auto Chamber=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Chamber.stl");
    //auto MgF2Lens=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/MgF2Lens.stl");
    //auto MgF2Lens=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/RotatedLens.stl");
    auto MgF2Window=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/MgF2Window.stl");
    auto AnodeTube=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/AnodeTube.stl");
    auto CathodeTube=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/CathodeTube.stl");
    //auto Peeks=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Peeks.stl");
    //auto Brackets=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Brackets.stl");
    auto Pmt=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/PMT.stl");
    auto Camera=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Camera.stl");
    auto Anode_Vacuum=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Anode_Vacuum.stl");
    auto Cathode_Vacuum=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Cathode_Vacuum.stl");
    auto Gas_Lens=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Gas_Lens.stl");
    auto Gas_Window=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Gas_Window.stl");

    //auto FieldCage_solid=FieldCage->GetSolid();
    //auto Meshes_solid=Meshes->GetSolid();
    auto Needle4_solid=Needle4->GetTessellatedSolid();
    auto Needle9_solid=Needle9->GetTessellatedSolid();
    auto Needle14_solid=Needle14->GetTessellatedSolid();
    auto Chamber_solid=Chamber->GetSolid();
    //auto MgF2Lens_solid=MgF2Lens->GetSolid();
    auto MgF2Window_solid=MgF2Window->GetSolid();
    auto AnodeTube_solid=AnodeTube->GetSolid();
    auto CathodeTube_solid=CathodeTube->GetSolid();
    //auto Peeks_solid=Peeks->GetSolid();
    //auto Brackets_solid=Brackets->GetSolid();
    auto Pmt_solid=Pmt->GetSolid();
    auto Camera_solid=Camera->GetSolid();
    auto AnodeVacuum_solid=Anode_Vacuum->GetSolid();
    auto CathodeVacuum_solid=Cathode_Vacuum->GetSolid();
    auto Gas_Lens_solid=Gas_Lens->GetSolid();
    auto Gas_Window_solid=Gas_Window->GetSolid();



    // Create Logical Space
    auto Needle4_logic= new G4LogicalVolume(Needle4_solid,Steel,"Needle4cm_logic");
    auto Needle9_logic=new G4LogicalVolume(Needle9_solid,Steel,"Needle9cm_logic");
    auto Needle14_logic=new G4LogicalVolume(Needle14_solid,Steel,"Needle14cm_logic");
    auto Chamber_logic=new G4LogicalVolume(Chamber_solid,Steel,"Chamber_logic");
    //auto FieldCage_logic=new G4LogicalVolume(FieldCage_solid,Steel,"FieldCage_logic");
    //auto MgF2Lens_logic=new G4LogicalVolume(MgF2Lens_solid,MgF2,"MgF2Lens_logic");
    auto AnodeTube_logic=new G4LogicalVolume(AnodeTube_solid,Steel,"AnodeTube_logic");
    //auto Peeks_logic=new G4LogicalVolume(Peeks_solid,PEEK,"Peeks_logic");
    //auto Brackets_logic=new G4LogicalVolume(Brackets_solid,materials::HDPE(),"Brackets_logic");
    //auto Meshes_logic=new G4LogicalVolume(Meshes_solid,Steel,"Mesh_logic");
    auto MgF2Window_logic=new G4LogicalVolume(MgF2Window_solid,gxe,"MgF2Window_logic");


    auto CathodeTube_logic=new G4LogicalVolume(CathodeTube_solid,Steel,"CathodeTube_logic");
    auto Pmt_logic=new G4LogicalVolume(Pmt_solid,MgF2,"Pmt_logic");
    auto Camera_logic=new G4LogicalVolume(Camera_solid,MgF2,"Camera_logic");
    auto AnodeVacuum_logic=new G4LogicalVolume(AnodeVacuum_solid,vacuum,"AnodeVacuum_logic");
    auto Gas_Lens_logic=new G4LogicalVolume(Gas_Lens_solid,gxe,"GasLens_logic");
    auto Gas_Window_logic=new G4LogicalVolume(Gas_Window_solid,gxe,"GasWindow_logic");
    auto CathodeVacuum_logic=new G4LogicalVolume(CathodeVacuum_solid,vacuum,"CathodeVacuum_logic");
    G4double Offset=0.2*cm;
    // Placement of the Items
    auto labPhysical = new G4PVPlacement(0, G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0, false,0,false);

    auto Chamber_physical=new G4PVPlacement(0,G4ThreeVector(),Chamber_logic,Chamber_solid->GetName(),lab_logic_volume,0,0,false);
    auto gas_pyhsical=new G4PVPlacement(0,G4ThreeVector(0,0,0),gas_logic,gas_logic->GetName(),lab_logic_volume,0,0,false);

    //auto FieldCage_physical=new G4PVPlacement(0,G4ThreeVector(),FieldCage_logic,FieldCage_solid->GetName(),gas_logic,false,0,false);

    // Peeks and Brackets
    //auto Peeks_physical=new G4PVPlacement(0,G4ThreeVector(),Peeks_logic,Peeks_solid->GetName(),gas_logic,0,0,false);
    //auto Brackets_physical=new G4PVPlacement(0,G4ThreeVector(),Brackets_logic,Brackets_solid->GetName(),gas_logic,0,0,false);
    // Gas Filling the Gaps
    auto Gas_Lens_pysical=new G4PVPlacement(0,G4ThreeVector(0,0,-1*mm/2-Offset+1*mm/2),Gas_Lens_logic,Gas_Lens_solid->GetName(),lab_logic_volume,0,0,false);
    //auto Gas_Lens_pysical=new G4PVPlacement(0,G4ThreeVector(0,0,1.5*mm/2+-Offset),Gas_Lens_logic,Gas_Lens_solid->GetName(),lab_logic_volume,0,0,false);
    auto Gas_Window_pysical=new G4PVPlacement(0,G4ThreeVector(0,0,1.5*mm/2),Gas_Window_logic,Gas_Window_solid->GetName(),lab_logic_volume,0,0,false);

    //Lens and Window
    //auto MgF2Lens_physical=new G4PVPlacement(0,G4ThreeVector(0,0,-Offset-0.1*cm),MgF2Lens_logic,MgF2Lens_solid->GetName(),lab_logic_volume,0,0,false);
    auto MgF2Window_physical=new G4PVPlacement(0,G4ThreeVector(),MgF2Window_logic,MgF2Window_solid->GetName(),lab_logic_volume,0,0,false);
    G4VPhysicalVolume* lensPhysical = new G4PVPlacement(0, G4ThreeVector(-1.1*mm,0.2*mm , 23.03*cm+ maxLensLength/2.0-Offset-2*mm/2), lensLogical,"MgF2_WINDOW1", lab_logic_volume,false, 0, false);

    // Meshes
    //auto Mesh_pysical=new G4PVPlacement(0,G4ThreeVector(),Meshes_logic,Meshes_solid->GetName(),gas_logic,0,0,false);

    // Image Plane
    auto Pmt_physical=new G4PVPlacement(0,G4ThreeVector(),Pmt_logic,Pmt_solid->GetName(),AnodeVacuum_logic,0,0,false);
    auto Camera_physical=new G4PVPlacement(0,G4ThreeVector(0,0,31*cm-Offset/2+1*mm/2),camLogical,camLogical->GetName(),CathodeVacuum_logic,0,0,false);
    //auto Camera_physical=new G4PVPlacement(0,G4ThreeVector(0,0,1.1*cm),Camera_logic,Camera_solid->GetName(),CathodeVacuum_logic,0,0,false);

    // Steel Tubes
    auto AnodeTube_physical=new G4PVPlacement(0,G4ThreeVector(),AnodeTube_logic,AnodeTube_solid->GetName(),AnodeVacuum_logic,0,0,false);
    auto CathodeTube_physical=new G4PVPlacement(0,G4ThreeVector(0,0,-Offset-0.1*cm),CathodeTube_logic,CathodeTube_solid->GetName(),CathodeVacuum_logic,0,0,false);
    // Vacuum in the Tubes
    auto AnodeVacuum_physical=new G4PVPlacement(0,G4ThreeVector(0,0,1*mm),AnodeVacuum_logic,AnodeVacuum_solid->GetName(),lab_logic_volume,0,0,false);
    auto CathodeVacuum_physical=new G4PVPlacement(0,G4ThreeVector(0,0,-Offset-0.1*cm),CathodeVacuum_logic,CathodeVacuum_solid->GetName(),lab_logic_volume,0,0,false);

    // Needle Placement
    auto Needle4_physical= new G4PVPlacement(0,G4ThreeVector(),Needle4_logic,Needle4_solid->GetName(),gas_logic,0,0,false);
    auto Needle9_physical= new G4PVPlacement(0,G4ThreeVector(),Needle9_logic,Needle9_solid->GetName(),gas_logic,0,0,false);
    auto Needle14_physical=new G4PVPlacement(0,G4ThreeVector(),Needle14_logic,Needle14_solid->GetName(),gas_logic,0,0,false);



    //////////////////////////////////////////

    // G4double EL_pos=-FielCageGap-ElGap_;
    G4double EL_pos=-10.98*cm;
    FielCageGap=21.26*cm;


    // These are mainly for illustration
    // EL Region
    G4Tubs* EL_solid = new G4Tubs("EL_GAP_", 0., Active_diam/2.,ElGap_/2 , 0., twopi);
    G4LogicalVolume* EL_logic = new G4LogicalVolume(EL_solid, gxe, "EL_GAP_");

    // EL_Gap
    new G4PVPlacement(0, G4ThreeVector(0.,0.,EL_pos-fOffset/2),EL_logic,EL_solid->GetName(),gas_logic, 0,0, false);
    // FieldCage -- needs to be updated to rings and PEEK rods
    G4Tubs* FieldCage_Solid =new G4Tubs("FIELDCAGE_", 0., Active_diam/2.,FielCageGap/2 , 0., twopi);
    G4LogicalVolume* FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "FIELDCAGE_");
    G4VPhysicalVolume * FieldCage_Phys=new G4PVPlacement(0,G4ThreeVector(0,0,-fOffset/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,false);


    /// Reflections from steel
    G4OpticalSurface * OpSteelSurf=new G4OpticalSurface("SteelSurface",unified,polished,dielectric_metal);
    OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());
    new G4LogicalBorderSurface("SteelSurface_CathodeChamber",CathodeVacuum_physical,Chamber_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_AnodeChamber",AnodeVacuum_physical,Chamber_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Chamber",gas_pyhsical,Chamber_physical,OpSteelSurf);
    //new G4LogicalBorderSurface("SteelSurface_Rings",gas_pyhsical,FieldCage_physical,OpSteelSurf);
    //new G4LogicalBorderSurface("SteelSurface_RingsAndMesh",gas_pyhsical,Mesh_pysical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Anode",AnodeVacuum_physical,AnodeTube_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Cathode",CathodeVacuum_physical,CathodeTube_physical,OpSteelSurf);


    // Reflection from needle
    new G4LogicalBorderSurface("SteelSurface_Needle4",gas_pyhsical,Needle4_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Needle9",gas_pyhsical,Needle9_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Needle14",gas_pyhsical,Needle14_physical,OpSteelSurf);

    G4OpticalSurface* opXenon_Glass2 = new G4OpticalSurface("XenonLensSurface");
    opXenon_Glass2->SetModel(glisur);                  // SetModel
    opXenon_Glass2->SetType(dielectric_dielectric);   // SetType
    opXenon_Glass2->SetFinish(polished);                 // SetFinish
    opXenon_Glass2->SetPolish(0.0);

    new G4LogicalBorderSurface("XenonLensSurface",Gas_Lens_pysical,lensPhysical,opXenon_Glass2);
    new G4LogicalBorderSurface("XenonLensSurface",Gas_Window_pysical,MgF2Window_physical,opXenon_Glass2);

    // For Camera and PMT assuming perfect detection
    G4OpticalSurface * CamSurf=new G4OpticalSurface("CamSurfaces",unified,ground,dielectric_metal);
    CamSurf->SetPolish(0);

    new G4LogicalBorderSurface("CameraSurface",CathodeVacuum_physical,Camera_physical,CamSurf);
    new G4LogicalBorderSurface("CameraSurface",AnodeVacuum_physical,Pmt_physical,CamSurf);


    // Call the Needles
    if(!HideSourceHolder_){
        // Particle Source Holder

        //Grab Points needed to sample from each or all the needles
        Sampler->SampleFromFacet(Needle4_solid);
        Sampler->SampleFromFacet(Needle9_solid);
        Sampler->SampleFromFacet(Needle14_solid);

       // This takes account of any shifting or rotation happens

        Sampler->SaveAllPointsToOneFile();
        G4VisAttributes *needlevis=new G4VisAttributes(G4Colour(1,1,1));
        needlevis->SetForceSolid(true);
        Needle9_logic->SetVisAttributes(needlevis);
        Needle4_logic->SetVisAttributes(needlevis);
        Needle14_logic->SetVisAttributes(needlevis);
    }

    // Visuals
    AssignVisuals();

    //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
    G4Region* regionGas = new G4Region("GasRegion");
    regionGas->AddRootLogicalVolume(gas_logic);
    OldCRAB *geo=new OldCRAB(fGasModelParameters);
    geo->SetMotherLab(lab_logic_volume);
    geo->SetOffset(0.2*cm);
    //geo->Construct();

#ifdef With_Opticks
    std::cout <<"Setting our detector geometry with opticks" <<std::endl;
    cudaDeviceReset();
    G4CXOpticks::SetGeometry(labPhysical);

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


void DetectorConstruction::AssignVisuals() {
    // Chamber
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
    // Lab
    G4LogicalVolume* Lab = lvStore->GetVolume("LAB");
    G4VisAttributes *LabVa=new G4VisAttributes(G4Colour(2,2,2));
    LabVa->SetForceWireframe(false);
    //Chamber
    G4LogicalVolume* Chamber = lvStore->GetVolume("Chamber_logic");
    G4VisAttributes *ChamberVa=new G4VisAttributes(G4Colour(1,1,1));
    ChamberVa->SetForceCloud(true);
    //ChamberVa->SetForceSolid(true);
    Chamber->SetVisAttributes(ChamberVa);


    //GAS
    G4LogicalVolume* Gas = lvStore->GetVolume("GAS_");
    G4VisAttributes *GasVa=new G4VisAttributes(colours::YellowAlpha());
    GasVa->SetForceCloud(true);
    Gas->SetVisAttributes(GasVa);

    G4LogicalVolume* Gas_Lens = lvStore->GetVolume("GasLens_logic");
    G4LogicalVolume* Gas_Window = lvStore->GetVolume("GasWindow_logic");
       Gas_Lens->SetVisAttributes(GasVa);
    Gas_Window->SetVisAttributes(GasVa);
    //PMT TUBE AND PMT BLOCK


    //Vacuum
    G4LogicalVolume* AnodeTubeVacuum = lvStore->GetVolume("AnodeVacuum_logic");
    G4LogicalVolume* CathodeTubeVacuum = lvStore->GetVolume("CathodeVacuum_logic");
    G4VisAttributes *VacumVis=new G4VisAttributes(colours::LillaAlpha());
    VacumVis->SetForceCloud(true);
    AnodeTubeVacuum->SetVisAttributes(VacumVis);
    CathodeTubeVacuum->SetVisAttributes(VacumVis);


    // Any Steel in Field Cage
    /*G4LogicalVolume* FRLog = lvStore->GetVolume("FieldCage_logic");
    G4VisAttributes FReVis= G4VisAttributes(colours::WhiteAlpha());
    FReVis.SetForceSolid(true);
    FRLog->SetVisAttributes(FReVis);

    // Brackets
    G4LogicalVolume* BracketLog = lvStore->GetVolume("Brackets_logic");
    G4VisAttributes BracketVis=colours::DirtyWhite();
    BracketVis.SetForceSolid(true);
    BracketLog->SetVisAttributes(BracketVis);


    // PEEK
    G4LogicalVolume* PEEKLog = lvStore->GetVolume("Peeks_logic");
    G4VisAttributes PEEKVis=colours::YellowAlpha();
    PEEKVis.SetForceSolid(true);
    PEEKLog->SetVisAttributes(PEEKVis);

 */

    //PMT TUBE AND PMT BLOCK
    G4LogicalVolume * AnodeTube=lvStore->GetVolume("AnodeTube_logic");
    G4LogicalVolume * CathodeTube=lvStore->GetVolume("CathodeTube_logic");
    AnodeTube->SetVisAttributes(ChamberVa);
    CathodeTube->SetVisAttributes(ChamberVa);

    //G4LogicalVolume * Meshes=lvStore->GetVolume("Mesh_logic");
    //Meshes->SetVisAttributes(FReVis);
   // PmttubeVacuumVis.SetForceCloud(true);




    //MgF2Window
    //G4LogicalVolume* lensLogical = lvStore->GetVolume("MgF2Lens_logic");
    G4LogicalVolume* WindowLogic = lvStore->GetVolume("MgF2Window_logic");
    G4VisAttributes  MgF2LensVis=colours::DarkGreen();
    //MgF2LensVis.SetForceSolid(true);
    //lensLogical->SetVisAttributes(MgF2LensVis);
    WindowLogic->SetVisAttributes(MgF2LensVis);


    //MgF2Window
    G4LogicalVolume* lensLogical2 = lvStore->GetVolume("Lens2");
    G4VisAttributes  MgF2LensVis2=colours::Red();
    MgF2LensVis2.SetForceSolid(true);
    lensLogical2->SetVisAttributes(MgF2LensVis2);
    // Camera
    G4LogicalVolume* CAMLog = lvStore->GetVolume("Camera_logic");
    G4LogicalVolume* PMT = lvStore->GetVolume("Pmt_logic");
    G4VisAttributes CAMVis=colours::DarkGreen();
    CAMVis.SetForceSolid(true);
    CAMLog->SetVisAttributes(CAMVis);
    PMT->SetVisAttributes(CAMVis);

    // EL-Region
    G4LogicalVolume * ELLogic=lvStore->GetVolume("EL_GAP_");
    G4VisAttributes ELVis=colours::BlueAlpha();
    ELVis.SetForceCloud(true);
    ELLogic->SetVisAttributes(ELVis);

    // FieldCage
    G4LogicalVolume * FieldCage=lvStore->GetVolume("FIELDCAGE_");
    G4VisAttributes FielCageVis=colours::Red();
    FielCageVis.SetForceCloud(true);
    FieldCage->SetVisAttributes(FielCageVis);

    Lab->SetVisAttributes(G4VisAttributes::GetInvisible());


    }
