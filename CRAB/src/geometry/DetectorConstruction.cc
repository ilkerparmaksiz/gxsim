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
    chamber_diam   (16.4 * cm),
    chamber_length (43.18 * cm), // Config files vary
    chamber_thickn (7. * mm),
    SourceEn_offset (5.7 *cm),
    SourceEn_diam   (1.0 * cm),
    SourceEn_length (1 * cm),
    SourceEn_thickn (2. * mm),
    SourceEn_holedia (5. * mm),
    gas_pressure_(10 * bar),
    vtx_(0*cm,-1.6*cm,-5.25*cm),
    Active_diam(8.6 * cm),
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
    G4LogicalVolume * lab_logic_volume = new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;


    // Detector Gas Region
    //
    G4Tubs* gas_solid =new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2. + chamber_thickn, 0., twopi);
    gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

    // Optical Properties Assigned here
    MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());
    vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
    gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68,sc_yield_,e_lifetime_));

    //std::cout<<filename<<std::endl;

    // Import Detector Geometry from an STL

    auto start = std::chrono::high_resolution_clock::now();
    auto FieldCage=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/FieldCage_Steel.stl");
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout<<"Time it TOOK to read this file is --  > " <<duration.count()<<std::endl;
    start = std::chrono::high_resolution_clock::now();

    auto Needle4=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Needle4cm.stl");
    auto Needle9=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Needle9cm.stl");
    auto Needle14=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Needle14cm.stl");
    auto Chamber=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Chamber.stl");
    auto LensAndWindow=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/MGF2.stl");
    auto SteelTubes=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/SteelTubes.stl");
    auto Peeks=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Peeks.stl");
    auto Brackets=CADMesh::TessellatedMesh::FromSTL(crabpath+"data/Brakets.stl");

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    Chamber->SetOffset(G4ThreeVector(x, y, z));
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout<<"Time it TOOK to read this file is --  > " <<duration.count()<<std::endl;
    // Collect the Solids
    start = std::chrono::high_resolution_clock::now();
    //auto FieldCage_solid=FieldCage->G();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout<<"Time it TOOK to GeT this SOLID is --  > " <<duration.count()<<std::endl;
    start = std::chrono::high_resolution_clock::now();

    auto Needle4_solid=Needle4->GetTessellatedSolid();
    auto Needle9_solid=Needle9->GetTessellatedSolid();
    auto Needle14_solid=Needle14->GetTessellatedSolid();
    auto Chamber_solid=Chamber->GetTessellatedSolid();
    auto LensAndWindow_solid=LensAndWindow->GetTessellatedSolid();
    auto SteelTubes_solid=SteelTubes->GetTessellatedSolid();
    auto Peeks_solid=Peeks->GetTessellatedSolid();
    auto Brackets_solid=Brackets->GetTessellatedSolid();
    auto FieldCage_solid=FieldCage->GetTessellatedSolid();

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout<<"Time it TOOK to GeT this SOLIDs is --  > " <<duration.count()<<std::endl;
    // Create Logical Space
    auto Needle4_logic= new G4LogicalVolume(Needle4_solid,Steel,"Needle4cm_logic");
    auto Needle9_logic=new G4LogicalVolume(Needle9_solid,Steel,"Needle9cm_logic");
    auto Needle14_logic=new G4LogicalVolume(Needle14_solid,Steel,"Needle14cm_logic");
    auto Chamber_logic=new G4LogicalVolume(Chamber_solid,Steel,"Chamber_logic");
    auto FieldCage_logic=new G4LogicalVolume(FieldCage_solid,Steel,"FieldCage_logic");
    auto LensAndWindow_logic=new G4LogicalVolume(LensAndWindow_solid,MgF2,"LensAndWindow_logic");
    auto SteelTubes_logic=new G4LogicalVolume(SteelTubes_solid,Steel,"SteelTubes_logic");
    auto Peeks_logic=new G4LogicalVolume(Peeks_solid,PEEK,"Peeks_logic");
    auto Brackets_logic=new G4LogicalVolume(Brackets_solid,PEEK,"Brackets_logic");

    // Placement of the Items

    auto Chamber_physical=new G4PVPlacement(0,G4ThreeVector(),Chamber_logic,Chamber_solid->GetName(),lab_logic_volume,0,0,false);
    auto FieldCage_physical=new G4PVPlacement(0,G4ThreeVector(),FieldCage_logic,FieldCage_solid->GetName(),Chamber_logic,0,0,false);
    auto LensAndWindow_physical=new G4PVPlacement(0,G4ThreeVector(),LensAndWindow_logic,LensAndWindow_solid->GetName(),Chamber_logic,0,0,false);
    auto SteelTubes_physical=new G4PVPlacement(0,G4ThreeVector(),SteelTubes_logic,SteelTubes_solid->GetName(),lab_logic_volume,0,0,false);
    auto Peeks_physical=new G4PVPlacement(0,G4ThreeVector(),Peeks_logic,Peeks_solid->GetName(),Chamber_logic,0,0,false);
    auto Brackets_physical=new G4PVPlacement(0,G4ThreeVector(),Brackets_logic,Brackets_solid->GetName(),Chamber_logic,0,0,false);



    //////////////////////////////////////////
    G4double FieldCagePos=-0*cm;
    // G4double EL_pos=-FielCageGap-ElGap_;
    G4double EL_pos=-10.98*cm;
    // FielCageGap=(16.03+2.955)*cm;
    FielCageGap=21.26*cm;


    // EL Region
    G4Tubs* EL_solid = new G4Tubs("EL_GAP", 0., Active_diam/2.,ElGap_/2 , 0., twopi);
    G4LogicalVolume* EL_logic = new G4LogicalVolume(EL_solid, gxe, "EL_GAP");

    // EL_Gap
    new G4PVPlacement(0, G4ThreeVector(0.,0.,EL_pos),EL_logic,EL_solid->GetName(),gas_logic, 0,0, false);




    // FieldCage -- needs to be updated to rings and PEEK rods
    G4Tubs* FieldCage_Solid =new G4Tubs("FIELDCAGE", 0., Active_diam/2.,FielCageGap/2 , 0., twopi);
    G4LogicalVolume* FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "FIELDCAGE");

    G4VPhysicalVolume * FieldCage_Phys=new G4PVPlacement(0,G4ThreeVector(0,0,FieldCagePos/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,false);

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
    G4LogicalVolume * pmt1_logic=pmt1_->GetLogicalVolume();
    G4LogicalVolume * pmt2_logic=pmt2_->GetLogicalVolume();


    // PMT1 and PMT3

    // // PMTs
    G4double PMT_offset=0.2*cm;
    G4double PMT_pos=(chamber_length/2)+chamber_thickn+(pmt1_->Length()/2)+MgF2_window_thickness_+PMT_offset;
    G4RotationMatrix* pmt1rotate = new G4RotationMatrix();
    pmt1rotate->rotateX(180.*deg);
    // pmt1rotate->rotateX(90.*deg);



    // CAMERA WINDOW
    G4double camHalfLength=0.5*mm;
    G4double camRadius= 12.7*mm;
    G4VSolid* camSolid = new G4Tubs("camWindow",0.,camRadius,camHalfLength,0.,twopi);
    G4LogicalVolume* camLogical = new G4LogicalVolume(camSolid,MgF2,"camLogical");



    auto labPhysical = new G4PVPlacement(0, G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0, false,0, false);


    // Xenon Gas in Active Area and Non-Active Area
    G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_logic, gas_solid->GetName(),lab_logic_volume, false, 0, false);

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




    /// OpticalSurface
    G4OpticalSurface * OpSteelSurf=new G4OpticalSurface("SteelSurface",unified,polished,dielectric_metal);
    OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());
    new G4LogicalBorderSurface("SteelSurface_Chamber",labPhysical,Chamber_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_RingsAndMesh",labPhysical,FieldCage_physical,OpSteelSurf);
    new G4LogicalBorderSurface("SteelSurface_ImagingDeviceHolders",labPhysical,FieldCage_physical,OpSteelSurf);


    // Camera
    G4OpticalSurface* opXenon_Glass = new G4OpticalSurface("XenonCamSurface");
    opXenon_Glass->SetModel(unified);                  // SetModel
    opXenon_Glass->SetType(dielectric_metal);   // SetType
    opXenon_Glass->SetFinish(polished);                 // SetFinish
    //new G4LogicalBorderSurface("XenonCamSurface",PMT_Tube_Vacuum_Phys0, camPhysical,opXenon_Glass);


    // Lens
    G4OpticalSurface* opXenon_Glass2 = new G4OpticalSurface("XenonLensSurface");
    opXenon_Glass2->SetModel(glisur);                  // SetModel
    opXenon_Glass2->SetType(dielectric_dielectric);   // SetType
    opXenon_Glass2->SetFinish(polished);                 // SetFinish
    opXenon_Glass2->SetPolish(0.0);
    new G4LogicalBorderSurface("XenonLensSurface",labPhysical,LensAndWindow_physical,opXenon_Glass2);
    //new G4LogicalBorderSurface("XenonPMTsSurface",gas_phys,PMTWindow,opXenon_Glass2);


    // Call the Needles
    if(!HideSourceHolder_){
        // Particle Source Holder

        if(!HideCollimator_) {
            //new G4PVPlacement(NeedleRotate,CollPosition,Coll_Logic,CollimatorWithBlock->GetName(),gas_logic,true,0,false);
        }


        auto Needle4_physical= new G4PVPlacement(0,G4ThreeVector(),Needle4_logic,Needle4_solid->GetName(),FieldCage_logic,0,0,false);
        auto Needle9_physical= new G4PVPlacement(0,G4ThreeVector(),Needle9_logic,Needle9_solid->GetName(),FieldCage_logic,0,0,false);
        auto Needle14_physical=new G4PVPlacement(0,G4ThreeVector(),Needle14_logic,Needle14_solid->GetName(),FieldCage_logic,0,0,false);
        Sampler->FaceTransform(Needle4_physical);
        Sampler->FaceTransform(Needle9_physical);
        Sampler->FaceTransform(Needle14_physical);



        G4VisAttributes *needlevis=new G4VisAttributes(G4Colour(1,1,1));
        needlevis->SetForceSolid(true);
        Needle9_logic->SetVisAttributes(needlevis);
        Needle4_logic->SetVisAttributes(needlevis);
        Needle14_logic->SetVisAttributes(needlevis);

        new G4LogicalBorderSurface("SteelSurface_Needle4",labPhysical,Needle4_physical,OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_Needle9",labPhysical,Needle9_physical,OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_Needle14",labPhysical,Needle14_physical,OpSteelSurf);

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
    G4VisAttributes *GasVa=new G4VisAttributes(colours::WhiteAlpha());
    GasVa->SetForceCloud(true);
    Gas->SetVisAttributes(GasVa);

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
    G4VisAttributes BracketVis=G4VisAttributes(G4Colour(0,0,0));
    BracketVis.SetForceSolid(true);
    BracketLog->SetVisAttributes(BracketVis);


    // PEEK
    G4LogicalVolume* PEEKLog = lvStore->GetVolume("Peeks_logic");
    G4VisAttributes PEEKVis=colours::YellowAlpha();
    PEEKVis.SetForceSolid(true);
    PEEKLog->SetVisAttributes(PEEKVis);



    //PMT TUBE AND PMT BLOCK
    G4LogicalVolume * SteelTubes=lvStore->GetVolume("SteelTubes_logic");
    SteelTubes->SetVisAttributes(ChamberVa);
    /*G4LogicalVolume * PmttubeBlockLog0=lvStore->GetVolume("PMT_TUBE_BLOCK0");
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
    G4LogicalVolume* lensLogical = lvStore->GetVolume("LensAndWindow_logic");
    G4VisAttributes  MgF2LensVis=colours::DarkGreen();
    MgF2LensVis.SetForceSolid(true);
    lensLogical->SetVisAttributes(MgF2LensVis);


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
