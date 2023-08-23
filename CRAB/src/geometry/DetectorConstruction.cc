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
#include <chrono>

// Lets use tetgen
#ifdef tetgen
#define USE_CADMESH_TETGEN 1
#endif
#include "CADMesh.hh"

DetectorConstruction::DetectorConstruction(GasModelParameters* gmp) :
    fGasModelParameters(gmp),
    checkOverlaps(1),
    temperature(300*kelvin), // temperature
    Lab_size(3. *m),
    chamber_diam   (19.4 * cm),
    chamber_length (41.8* cm), // Config files vary
    //chamber_length (42.5* cm), // Config files vary
    chamber_thickn (7. * mm),
    SourceEn_offset (5.7 *cm),
    SourceEn_diam   (1.0 * cm),
    SourceEn_length (1 * cm),
    SourceEn_thickn (2. * mm),
    SourceEn_holedia (5. * mm),
    gas_pressure_(10 * bar),
    vtx_(0*cm,-1.6*cm,-5.25*cm),
    Active_diam(7.2 * cm),
    sc_yield_(25510./MeV),
    e_lifetime_(1000. * ms),
    pmt_hole_length_ (18.434 * cm),
    MgF2_window_thickness_ (6. * mm),
    MgF2_window_diam_ (16.5 * mm),
    wndw_ring_stand_out_ (1.5 * mm), //how much the ring around sapph windows stands out of them
    pedot_coating_thickness_ (200. * nanometer), // copied from NEW
    optical_pad_thickness_ (1. * mm), // copied from NEW
    pmt_base_diam_ (47. * mm),
    pmt_base_thickness_ (5. * mm),
    HideSourceHolder_(false),
    max_step_size_(1.*mm),
    ElGap_(7*mm),
    ELyield_(970/cm),
    PMT1_Pos_(2.32*cm),
    PMT3_Pos_(3.52*cm),
    HideCollimator_(true)
{
    detectorMessenger = new DetectorMessenger(this);
    Sampler=new util::SampleFromSurface("Needle");
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

    G4Tubs* gas_solid =new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2. + chamber_thickn+3*cm/2, 0., twopi);
    //G4Tubs* gas_solid =new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2. + chamber_thickn, 0., twopi);


    // Optical Properties Assigned here
    MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());
    vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
    gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68,sc_yield_,e_lifetime_));

    G4LogicalVolume * lab_logic_volume = new G4LogicalVolume(lab_solid_volume,air,lab_name) ; //G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR")

    gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");


    /// MgF2 window ///
    G4Tubs* MgF2_window_solid = new G4Tubs("MgF2_WINDOW", 0., MgF2_window_diam_/2.,
                                           (MgF2_window_thickness_ )/2., 0., twopi);
    G4LogicalVolume* MgF2_window_logic= new G4LogicalVolume(MgF2_window_solid, MgF2, "MgF2_WINDOW");

    // lens
    const G4double lensRcurve (2.83*cm); // radius of curvature of MgF2 Lens
    const G4ThreeVector posLensTubeIntersect (0.,0.,-lensRcurve);

    // Create lens from the intersection of a sphere and a cylinder
    G4double maxLensLength = 4*mm;
    G4Tubs* sLensTube = new G4Tubs("sLensSphereTube", 0, MgF2_window_diam_/2, maxLensLength, 0.,twopi); // 4 mm is the max lens length
    G4Orb* sLensOrb = new G4Orb("sLensSphere",lensRcurve);
    G4IntersectionSolid* sLens =  new G4IntersectionSolid("sLens",sLensTube,sLensOrb, 0, posLensTubeIntersect);

    // Lens logical
    G4LogicalVolume* lensLogical = new G4LogicalVolume(sLens, MgF2, "Lens");

    // Import Detector Geometry from an STL

    auto FieldCage=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/FieldRings.stl");
    //auto Meshes=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Meshes.stl");
    auto Needle4=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_4cm.stl");
    auto Needle9=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_9cm.stl");
    auto Needle14=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Needle_14cm.stl");
    auto Chamber=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Chamber.stl");
    auto MgF2Lens=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/MgF2Lens.stl");
    auto MgF2Window=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/MgF2Window.stl");
    auto AnodeTube=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/AnodeTube.stl");
    auto CathodeTube=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/CathodeTube.stl");
    //auto Peeks=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Peeks.stl");
    auto Brackets=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Brackets.stl");
    auto Pmt=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/PMT.stl");
    auto Camera=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Camera.stl");
    auto Anode_Vacuum=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Anode_Vacuum.stl");
    auto Cathode_Vacuum=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Cathode_Vacuum.stl");
    auto Gas_Lens=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Gas_Lens.stl");
    auto Gas_Window=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/CRAB_STL/Gas_Window.stl");

    auto FieldCage_solid=FieldCage->GetTessellatedSolid();
    //auto Meshes_solid=Meshes->GetTessellatedSolid();
    auto Needle4_solid=Needle4->GetTessellatedSolid();
    auto Needle9_solid=Needle9->GetTessellatedSolid();
    auto Needle14_solid=Needle14->GetTessellatedSolid();
    auto Chamber_solid=Chamber->GetTessellatedSolid();
    auto MgF2Lens_solid=MgF2Lens->GetTessellatedSolid();
    auto MgF2Window_solid=MgF2Window->GetTessellatedSolid();
    auto AnodeTube_solid=AnodeTube->GetTessellatedSolid();
    auto CathodeTube_solid=CathodeTube->GetTessellatedSolid();
    //auto Peeks_solid=Peeks->GetTessellatedSolid();
    auto Brackets_solid=Brackets->GetTessellatedSolid();
    auto Pmt_solid=Pmt->GetTessellatedSolid();
    auto Camera_solid=Camera->GetTessellatedSolid();
    auto AnodeVacuum_solid=Anode_Vacuum->GetTessellatedSolid();
    auto CathodeVacuum_solid=Cathode_Vacuum->GetTessellatedSolid();
    auto Gas_Lens_solid=Gas_Lens->GetTessellatedSolid();
    auto Gas_Window_solid=Gas_Window->GetTessellatedSolid();



    // Create Logical Space
    auto Needle4_logic= new G4LogicalVolume(Needle4_solid,Steel,"Needle4cm_logic");
    auto Needle9_logic=new G4LogicalVolume(Needle9_solid,Steel,"Needle9cm_logic");
    auto Needle14_logic=new G4LogicalVolume(Needle14_solid,Steel,"Needle14cm_logic");
    auto Chamber_logic=new G4LogicalVolume(Chamber_solid,Steel,"Chamber_logic");
    auto FieldCage_logic=new G4LogicalVolume(FieldCage_solid,Steel,"FieldCage_logic");
    auto MgF2Lens_logic=new G4LogicalVolume(MgF2Lens_solid,MgF2,"MgF2Lens_logic");
    auto AnodeTube_logic=new G4LogicalVolume(AnodeTube_solid,Steel,"AnodeTube_logic");
    //auto Peeks_logic=new G4LogicalVolume(Peeks_solid,PEEK,"Peeks_logic");
    auto Brackets_logic=new G4LogicalVolume(Brackets_solid,materials::HDPE(),"Brackets_logic");
    //auto Meshes_logic=new G4LogicalVolume(Meshes_solid,Steel,"Mesh_logic");
    auto MgF2Window_logic=new G4LogicalVolume(MgF2Window_solid,gxe,"MgF2Window_logic");


    auto CathodeTube_logic=new G4LogicalVolume(CathodeTube_solid,Steel,"CathodeTube_logic");
    auto Pmt_logic=new G4LogicalVolume(Pmt_solid,MgF2,"Pmt_logic");
    auto Camera_logic=new G4LogicalVolume(Camera_solid,MgF2,"Camera_logic");
    auto AnodeVacuum_logic=new G4LogicalVolume(AnodeVacuum_solid,vacuum,"AnodeVacuum_logic");
    auto Gas_Lens_logic=new G4LogicalVolume(Gas_Lens_solid,gxe,"GasLens_logic");
    auto Gas_Window_logic=new G4LogicalVolume(Gas_Window_solid,gxe,"GasWindow_logic");
    auto CathodeVacuum_logic=new G4LogicalVolume(CathodeVacuum_solid,vacuum,"CathodeVacuum_logic");

    // Placement of the Items
    auto labPhysical = new G4PVPlacement(0, G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0, false,0,false);

    auto Chamber_physical=new G4PVPlacement(0,G4ThreeVector(),Chamber_logic,Chamber_solid->GetName(),lab_logic_volume,0,0,false);
    //auto gas_pyhsical=new G4PVPlacement(0,G4ThreeVector(0,0,0.7*cm/2),gas_logic,gas_logic->GetName(),lab_logic_volume,0,0,false);
    auto gas_pyhsical=new G4PVPlacement(0,G4ThreeVector(0,0,-0.2*cm/2),gas_logic,gas_logic->GetName(),lab_logic_volume,0,0,false);

    auto FieldCage_physical=new G4PVPlacement(0,G4ThreeVector(),FieldCage_logic,FieldCage_solid->GetName(),gas_logic,false,0,false);

    // Peeks and Brackets
    //auto Peeks_physical=new G4PVPlacement(0,G4ThreeVector(),Peeks_logic,Peeks_solid->GetName(),gas_logic,0,0,false);
    auto Brackets_physical=new G4PVPlacement(0,G4ThreeVector(),Brackets_logic,Brackets_solid->GetName(),gas_logic,0,0,false);
    // Gas Filling the Gaps
    //auto Gas_Lens_pysical=new G4PVPlacement(0,G4ThreeVector(0,0,-1*mm/2),Gas_Lens_logic,Gas_Lens_solid->GetName(),lab_logic_volume,0,0,false);
    //auto Gas_Window_pysical=new G4PVPlacement(0,G4ThreeVector(0,0,1.5*mm/2),Gas_Window_logic,Gas_Window_solid->GetName(),lab_logic_volume,0,0,false);

    //Lens and Window
    auto MgF2Lens_physical=new G4PVPlacement(0,G4ThreeVector(),MgF2Lens_logic,MgF2Lens_solid->GetName(),lab_logic_volume,0,0,false);
    auto MgF2Window_physical=new G4PVPlacement(0,G4ThreeVector(),MgF2Window_logic,MgF2Window_solid->GetName(),lab_logic_volume,0,0,false);

    // Meshes
    //auto Meshes_physical=new G4PVPlacement(0,G4ThreeVector(),Meshes_logic,Meshes_solid->GetName(),gas_logic,0,0,false);

    // Image Plane
    auto Pmt_physical=new G4PVPlacement(0,G4ThreeVector(),Pmt_logic,Pmt_solid->GetName(),AnodeVacuum_logic,0,0,false);
    auto Camera_physical=new G4PVPlacement(0,G4ThreeVector(),Camera_logic,Camera_solid->GetName(),CathodeVacuum_logic,0,0,false);

    // Vacuum in the Tubes
    auto AnodeVacuum_physical=new G4PVPlacement(0,G4ThreeVector(0,0,1*mm),AnodeVacuum_logic,AnodeVacuum_solid->GetName(),lab_logic_volume,0,0,false);
    auto CathodeVacuum_physical=new G4PVPlacement(0,G4ThreeVector(),CathodeVacuum_logic,CathodeVacuum_solid->GetName(),lab_logic_volume,0,0,false);
    auto AnodeTube_physical=new G4PVPlacement(0,G4ThreeVector(),AnodeTube_logic,AnodeTube_solid->GetName(),lab_logic_volume,0,0,false);
    auto CathodeTube_physical=new G4PVPlacement(0,G4ThreeVector(),CathodeTube_logic,CathodeTube_solid->GetName(),lab_logic_volume,0,0,false);


    // Needle Placement
    auto Needle4_physical= new G4PVPlacement(0,G4ThreeVector(),Needle4_logic,Needle4_solid->GetName(),gas_logic,0,0,false);
    auto Needle9_physical= new G4PVPlacement(0,G4ThreeVector(),Needle9_logic,Needle9_solid->GetName(),gas_logic,0,0,false);
    auto Needle14_physical=new G4PVPlacement(0,G4ThreeVector(),Needle14_logic,Needle14_solid->GetName(),gas_logic,0,0,false);



    //// PMT Covering Tube ///
    // // PMTs
    G4double PMT_offset=0.2*cm;

    G4double PMT_pos=(chamber_length/2)+chamber_thickn+MgF2_window_thickness_+PMT_offset+2*cm;
    G4RotationMatrix* pmt1rotate = new G4RotationMatrix();
    pmt1rotate->rotateX(180.*deg);
    G4double offset=1.65*cm;
    //G4double PMT_Tube_Length1=MgF2_window_thickness_+(pmt1_->Length()+0.5*cm)/2 + offset -PMT_offset-0.05*cm ;
    G4double PMT_Tube_Length1=MgF2_window_thickness_+(0.5*cm)/2 + offset -PMT_offset-0.05*cm ;
    //G4double PMT_Tube_Length0=(17*cm+pmt1_->Length())/2 - 3.87*cm-PMT_offset;
    G4double PMT_Tube_Length0=(17*cm)/2 - 3.87*cm-PMT_offset;
    G4double PMT_Tube_Block_Thickness=0.2*cm;
    G4double LongPMTTubeOffset=7.5*cm-3.9*cm;
    G4double PMTTubeDiam=2.54*cm;

    /*
    // Vacuum for PMT TUBE0
    /// add 2.54*cm to config file
    G4Tubs * InsideThePMT_Tube_solid0=new G4Tubs("PMT_TUBE_VACUUM0",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Length0,0,twopi);
    G4LogicalVolume * InsideThePMT_Tube_Logic0=new G4LogicalVolume(InsideThePMT_Tube_solid0,vacuum,InsideThePMT_Tube_solid0->GetName());

    // Tube Close to EL
    G4Tubs * PMT_Tube_solid1=new G4Tubs("PMT_TUBE1",(PMTTubeDiam/2)+0.5*cm,(PMTTubeDiam/2)+0.7*cm,PMT_Tube_Length1,0,twopi);
    G4LogicalVolume * PMT_Tube_Logic1=new G4LogicalVolume(PMT_Tube_solid1,materials::Steel(),PMT_Tube_solid1->GetName());
    G4Tubs * PMT_Block_solid1=new G4Tubs("PMT_TUBE_BLOCK1",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Block_Thickness,0,twopi);
    G4LogicalVolume * PMT_Block_Logic=new G4LogicalVolume(PMT_Block_solid1,materials::Steel(),PMT_Block_solid1->GetName());

    // Vacuum for PMT TUBE1

    G4Tubs * InsideThePMT_Tube_solid1=new G4Tubs("PMT_TUBE_VACUUM1",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Length1,0,twopi);
    G4LogicalVolume * InsideThePMT_Tube_Logic1=new G4LogicalVolume(InsideThePMT_Tube_solid1,vacuum,InsideThePMT_Tube_solid1->GetName());

    // Tube Away from EL
    G4Tubs * PMT_Tube_solid0=new G4Tubs("PMT_TUBE0",(PMTTubeDiam/2)+0.5*cm,(PMTTubeDiam/2)+0.8*cm,PMT_Tube_Length0,0,twopi);
    G4LogicalVolume * PMT_Tube_Logic0=new G4LogicalVolume(PMT_Tube_solid0,materials::Steel(),PMT_Tube_solid0->GetName());

    G4Tubs * PMT_Block_solid0=new G4Tubs("PMT_TUBE_BLOCK0",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Block_Thickness,0,twopi);
    G4LogicalVolume * PMT_Block_Logic0=new G4LogicalVolume(PMT_Block_solid0,materials::Steel(),PMT_Block_solid0->GetName());

    // PMT Tubes
    G4VPhysicalVolume *PMT_Tube_Phys0=new G4PVPlacement(0,G4ThreeVector(0, 0, PMT_pos + LongPMTTubeOffset-10*mm/2),     PMT_Tube_Logic0,PMT_Tube_Logic0->GetName(),lab_logic_volume,false,0,false);
    G4VPhysicalVolume *PMT_Tube_Phys1=new G4PVPlacement(0,G4ThreeVector(0, 0, -(PMT_pos - PMT_offset) - offset),PMT_Tube_Logic1,PMT_Tube_Logic1->GetName(),lab_logic_volume,false,0,false);

    // PMT Tube Vacuum
    G4VPhysicalVolume *PMT_Tube_Vacuum_Phys0=new G4PVPlacement(0,G4ThreeVector(0, 0, PMT_pos+LongPMTTubeOffset-10*mm/2),   InsideThePMT_Tube_Logic0,"PMT_TUBE_VACUUM0",lab_logic_volume,false,0,false);
    G4VPhysicalVolume *PMT_Tube_Vacuum_Phys1=new G4PVPlacement(0,G4ThreeVector(0, 0, -(PMT_pos-PMT_offset)-offset),InsideThePMT_Tube_Logic1,"PMT_TUBE_VACUUM1",lab_logic_volume,false,0,false);

    // PMT Tube Block
    new G4PVPlacement(0,G4ThreeVector(0,0, PMT_pos-PMT_offset+PMT_Tube_Length0-PMT_Tube_Block_Thickness/2+LongPMTTubeOffset),PMT_Block_Logic0,PMT_Block_Logic0->GetName(),lab_logic_volume,false,0,false);
    new G4PVPlacement(0,G4ThreeVector(0,0, -(PMT_pos-PMT_offset+PMT_Tube_Length1-PMT_Tube_Block_Thickness/2)-offset),        PMT_Block_Logic,PMT_Block_Logic->GetName(),lab_logic_volume,false,1,false);
    */
    // PMTs
     //new G4PVPlacement(pmt1rotate,G4ThreeVector (0,0,((PMT3_Pos_)-pmt1_->Length()/2-PMT_Tube_Length1/2-MgF2_window_thickness_/2)),pmt1_logic,pmt1_->GetPMTName(),InsideThePMT_Tube_Logic0,true,0,false);
    //new G4PVPlacement(0,G4ThreeVector (0, 0., (PMT1_Pos_-pmt1_->Length()/2-MgF2_window_thickness_/2)),pmt2_logic,pmt2_->GetPMTName(),InsideThePMT_Tube_Logic1,true,0,false);
    G4double window_posz = chamber_length/2 + chamber_thickn;
    G4RotationMatrix * rt=new G4RotationMatrix ();
    rt->rotateY(180*degree);
    //G4VPhysicalVolume* lensPhysical = new G4PVPlacement(0, G4ThreeVector(-1*mm, -0.9*mm/2, window_posz+ maxLensLength/2.0+0.9*cm), lensLogical,"MgF2_LENS", lab_logic_volume,false, 0, false);
    //G4VPhysicalVolume* WindowPhysical = new G4PVPlacement(0, G4ThreeVector(-1*mm, -0.9*mm/2, -window_posz-1.75*cm), MgF2_window_logic,"MgF2_WINDOW", lab_logic_volume, false, 1, false);






    //////////////////////////////////////////
    G4double FieldCagePos=-0*cm;
    G4double Offset=-0.8*cm;
    // G4double EL_pos=-FielCageGap-ElGap_;
    G4double EL_pos=-10.98*cm;
    // FielCageGap=(16.03+2.955)*cm;
    FielCageGap=21.26*cm;


    // EL Region
    G4Tubs* EL_solid = new G4Tubs("EL_GAP", 0., Active_diam/2.,ElGap_/2 , 0., twopi);
    G4LogicalVolume* EL_logic = new G4LogicalVolume(EL_solid, gxe, "EL_GAP");

    // EL_Gap
    new G4PVPlacement(0, G4ThreeVector(0.,0.,EL_pos-Offset/2),EL_logic,EL_solid->GetName(),gas_logic, 0,0, true);




    // FieldCage -- needs to be updated to rings and PEEK rods
    G4Tubs* FieldCage_Solid =new G4Tubs("FIELDCAGE", 0., Active_diam/2.,FielCageGap/2 , 0., twopi);
    G4LogicalVolume* FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "FIELDCAGE");

    G4VPhysicalVolume * FieldCage_Phys=new G4PVPlacement(0,G4ThreeVector(0,0,FieldCagePos/2-Offset/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,false);

    // Radioactive Source Encloser
    // Source
    G4Tubs* SourceHolChamber_solid =new G4Tubs("SourceHolChamber", SourceEn_holedia/2, (SourceEn_diam/2. + SourceEn_thickn),(SourceEn_length/2. + SourceEn_thickn),0,twopi);
    G4LogicalVolume* SourceHolChamber_logic = new G4LogicalVolume(SourceHolChamber_solid,materials::Steel(), "SourceHolChamber_logic");

    G4Tubs* SourceHolChamberBlock_solid =new G4Tubs("SourceHolChBlock",0,(SourceEn_holedia/2),( SourceEn_thickn/2), 0.,twopi);
    G4LogicalVolume* SourceHolChamberBlock_logic = new G4LogicalVolume(SourceHolChamberBlock_solid,materials::Steel(), "SourceHolChBlock_logic");




    //Call Sampler for Needles
    Sampler->SampleFromFacet(Needle4_solid);
    Sampler->SampleFromFacet(Needle9_solid);
    Sampler->SampleFromFacet(Needle14_solid);
    auto SamplePoints= Sampler->getRawPoints();

    ///
    //Adding the PMTs in here
    pmt1_=new pmt::PmtR7378A();
    pmt2_=new pmt::PmtR7378A();
    pmt1_->SetPMTName("S2");
    pmt2_->SetPMTName("S1");
    pmt1_->Construct();
    pmt2_->Construct();


    // Adding Logical Volumes for PMTs
    //G4LogicalVolume * pmt1_logic=pmt1_->GetLogicalVolume();
    //G4LogicalVolume * pmt2_logic=pmt2_->GetLogicalVolume();




    // CAMERA WINDOW
    G4double camHalfLength=0.5*mm;
    G4double camRadius= 12.7*mm;
    G4VSolid* camSolid = new G4Tubs("camWindow",0.,camRadius,camHalfLength,0.,twopi);
    G4LogicalVolume* camLogical = new G4LogicalVolume(camSolid,MgF2,"camLogical");



    // Xenon Gas in Active Area and Non-Active Area
    //G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0.,Offset/6), gas_logic, gas_solid->GetName(),lab_logic_volume, false, 0, false);

    // FieldCage
    // G4double FieldCagePos=chamber_length/2-((129)*mm)-FielCageGap/2-ElGap_/2;
    // G4double EL_pos=chamber_length/2-FielCageGap/2-((326)*mm)-ElGap_/2;

    //G4double FieldCagePos=-0*cm;
    // G4double EL_pos=-FielCageGap-ElGap_;
    //G4double EL_pos=-10.98*cm;

    // G4VPhysicalVolume * FieldCage_Phys=new G4PVPlacement(0,G4ThreeVector(0,0,FieldCagePos/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,false)

    // PMTs
    // new G4PVPlacement(pmt1rotate,G4ThreeVector (0,0,((PMT3_Pos_)-pmt1_->Length()/2-PMT_Tube_Length1/2-MgF2_window_thickness_/2)),pmt1_logic,pmt1_->GetPMTName(),InsideThePMT_Tube_Logic0,true,0,false);
    //new G4PVPlacement(0,G4ThreeVector (0, 0., (PMT1_Pos_-pmt1_->Length()/2-MgF2_window_thickness_/2)),pmt2_logic,pmt2_->GetPMTName(),InsideThePMT_Tube_Logic1,true,0,false);


    // Place the camera Make camLogical mother and photocathode daughter
    // G4double LensFocalDist = 6.34*cm; // Got from trial and erro
    G4double ImageDist = 7.945*cm; // Got from trial and error
    //G4VPhysicalVolume* camPhysical= new G4PVPlacement(0,  G4ThreeVector (0,0, (chamber_length/2 + chamber_thickn + ImageDist) - PMT_pos-LongPMTTubeOffset),camLogical,"camPhysical",InsideThePMT_Tube_Logic0, false,0,false);




    /// Reflections from steel
    G4OpticalSurface * OpSteelSurf=new G4OpticalSurface("SteelSurface",unified,polished,dielectric_metal);
    OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());
    new G4LogicalBorderSurface("SteelSurface_Chamber",gas_pyhsical,Chamber_physical,OpSteelSurf);
    //new G4LogicalBorderSurface("SteelSurface_Chamber",labPhysical,Chamber_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_RingsAndMesh",gas_pyhsical,FieldCage_physical,OpSteelSurf);
    //new G4LogicalBorderSurface("SteelSurface_Anode",AnodeVacuum_physical,AnodeTube_physical,OpSteelSurf);
    //new G4LogicalBorderSurface("SteelSurface_Cathode",CathodeVacuum_physical,CathodeTube_physical,OpSteelSurf);

    // Reflection from needle
    new G4LogicalBorderSurface("SteelSurface_Needle4",gas_pyhsical,Needle4_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Needle9",gas_pyhsical,Needle9_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_Needle14",gas_pyhsical,Needle14_physical,OpSteelSurf);

    G4OpticalSurface* opXenon_Glass = new G4OpticalSurface("XenonSurface");
    opXenon_Glass->SetModel(glisur);                  // SetModel
    opXenon_Glass->SetType(dielectric_dielectric);   // SetType
    opXenon_Glass->SetFinish(polished);                 // SetFinish
    new G4LogicalBorderSurface("XenonSurfaceWindow",gas_pyhsical,MgF2Window_physical ,opXenon_Glass);
    new G4LogicalBorderSurface("XenonSurfaceLens",gas_pyhsical , MgF2Lens_physical, opXenon_Glass);
    //new G4LogicalBorderSurface("XenonSurfaceWindow",gas_pyhsical,lensPhysical ,opXenon_Glass);
    //new G4LogicalBorderSurface("XenonSurfaceLens",gas_pyhsical , WindowPhysical, opXenon_Glass);




    // Call the Needles
    if(!HideSourceHolder_){
        // Particle Source Holder

        if(!HideCollimator_) {
            //new G4PVPlacement(NeedleRotate,CollPosition,Coll_Logic,CollimatorWithBlock->GetName(),gas_logic,true,0,false);
        }

        Sampler->FaceTransform(Needle4_physical);
        Sampler->FaceTransform(Needle9_physical);
        Sampler->FaceTransform(Needle14_physical);



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
    Chamber->SetVisAttributes(ChamberVa);


    //GAS
    G4LogicalVolume* Gas = lvStore->GetVolume("GAS");
    G4VisAttributes *GasVa=new G4VisAttributes(colours::YellowAlpha());
    GasVa->SetForceCloud(true);
    Gas->SetVisAttributes(GasVa);

    G4LogicalVolume* Gas_Lens = lvStore->GetVolume("GasLens_logic");
    G4LogicalVolume* Gas_Window = lvStore->GetVolume("GasWindow_logic");
       Gas_Lens->SetVisAttributes(GasVa);
    Gas_Window->SetVisAttributes(GasVa);

    //Vacuum
    G4LogicalVolume* AnodeTubeVacuum = lvStore->GetVolume("AnodeVacuum_logic");
    G4LogicalVolume* CathodeTubeVacuum = lvStore->GetVolume("CathodeVacuum_logic");
    G4VisAttributes *VacumVis=new G4VisAttributes(colours::LillaAlpha());
    VacumVis->SetForceCloud(true);
    AnodeTubeVacuum->SetVisAttributes(VacumVis);
    CathodeTubeVacuum->SetVisAttributes(VacumVis);

    //Source Enclosure Related
    G4LogicalVolume* SourceHolder = lvStore->GetVolume("SourceHolChamber_logic");
    G4LogicalVolume* SourceHolderBlock = lvStore->GetVolume("SourceHolChBlock_logic");
    G4VisAttributes *SourceHolderVa=new G4VisAttributes(G4Colour(2,2,2));
    SourceHolderVa->SetForceSolid(true);


    // Any Steel in Field Cage
    G4LogicalVolume* FRLog = lvStore->GetVolume("FieldCage_logic");
    G4VisAttributes FReVis= G4VisAttributes(colours::WhiteAlpha());
    FReVis.SetForceSolid(true);
    FRLog->SetVisAttributes(FReVis);

    // Brackets
    G4LogicalVolume* BracketLog = lvStore->GetVolume("Brackets_logic");
    G4VisAttributes BracketVis=colours::DirtyWhite();
    BracketVis.SetForceSolid(true);
    BracketLog->SetVisAttributes(BracketVis);

    /*
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
    /*
    G4LogicalVolume * Meshes=lvStore->GetVolume("Mesh_logic");
    Meshes->SetVisAttributes(FReVis);
    G4LogicalVolume * PmttubeBlockLog0=lvStore->GetVolume("PMT_TUBE_BLOCK0");
    G4LogicalVolume * PmttubeLog1=lvStore->GetVolume("PMT_TUBE1");
    PmttubeLog1->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4LogicalVolume * PmttubeBlockLog1=lvStore->GetVolume("PMT_TUBE_BLOCK1");
    PmttubeBlockLog0->SetVisAttributes(ChamberVa);
    PmttubeBlockLog1->SetVisAttributes(ChamberVa);
    G4LogicalVolume * PmttubeVacuumLog1=lvStore->GetVolume("PMT_TUBE_VACUUM0");
    G4LogicalVolume * PmttubeVacuumLog2=lvStore->GetVolume("PMT_TUBE_VACUUM1");
    G4VisAttributes PmttubeVacuumVis=colours::DarkGreyAlpha();
    PmttubeVacuumVis.SetForceCloud(true);
    PmttubeVacuumLog1->SetVisAttributes(PmttubeVacuumVis);
    PmttubeVacuumLog2->SetVisAttributes(PmttubeVacuumVis);
    */

    //MgF2Window
    G4LogicalVolume* lensLogical = lvStore->GetVolume("MgF2Lens_logic");
    G4LogicalVolume* WindowLogic = lvStore->GetVolume("MgF2Window_logic");
    G4LogicalVolume* Window = lvStore->GetVolume("MgF2_WINDOW");
    G4LogicalVolume* Lens = lvStore->GetVolume("Lens");
    G4VisAttributes  MgF2LensVis=colours::DarkGreen();
    MgF2LensVis.SetForceSolid(true);
    lensLogical->SetVisAttributes(MgF2LensVis);
    WindowLogic->SetVisAttributes(MgF2LensVis);
    Lens->SetVisAttributes(MgF2LensVis);
    Window->SetVisAttributes(MgF2LensVis);


    // Camera
    G4LogicalVolume* CAMLog = lvStore->GetVolume("camLogical");
    G4VisAttributes CAMVis=colours::DarkRedAlpha();
    CAMVis.SetForceSolid(true);
    CAMLog->SetVisAttributes(CAMLog);

    // EL-Region
    G4LogicalVolume * ELLogic=lvStore->GetVolume("EL_GAP");
    G4VisAttributes ELVis=colours::BlueAlpha();
    ELVis.SetForceCloud(true);
    ELLogic->SetVisAttributes(ELVis);


    // FieldCage
    G4LogicalVolume * FieldCage=lvStore->GetVolume("FIELDCAGE");
    G4VisAttributes FielCageVis=colours::Red();
    FielCageVis.SetForceCloud(true);
    FieldCage->SetVisAttributes(FielCageVis);


    SourceHolder->SetVisAttributes(SourceHolderVa);
    SourceHolderBlock->SetVisAttributes(SourceHolderVa);
    Lab->SetVisAttributes(G4VisAttributes::GetInvisible());

    }
