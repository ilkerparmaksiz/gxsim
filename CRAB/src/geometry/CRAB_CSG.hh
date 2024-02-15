#ifndef OldCRAB_hh
#define OldCRAB_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorMessenger.hh"
#include "G4UserLimits.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4FieldManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UniformMagField.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4UnionSolid.hh"
#include "G4Region.hh"
#include "G4Orb.hh"
#include "GasModelParameters.hh"

#include "OpticalMaterialProperties.hh"
#include "MaterialsList.hh"
#include "PmtR7378A.hh"


class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UniformMagField;


using namespace std;
/*! \class  DetectorConstruction*/
/*! \brief class derived from G4VUserDetectorConstruction*/

class CRAB_CSG  {
 public:
    CRAB_CSG(GasModelParameters*);
    virtual ~CRAB_CSG();

    // Mandatory methods
    void  Construct();

    virtual void AssignVisuals();

    //Setters for the dimensions and environment variables of the setup
    inline void CheckOverlaps(G4bool co){checkOverlaps=co;};
    inline void SetGasPressure(G4double d){gas_pressure_=d;};
    inline void SetTemperature(G4double d){temperature=d;};
    inline void SetRegion(G4Region *region){rg=region;};
    inline void SetGasLogic(G4LogicalVolume * gas){gas_logic=gas;};

    inline G4double GetChamberR(){return chamber_diam/2.0/cm;};
    inline G4double GetChamberL(){return chamber_length/cm; }; 
    inline G4double GetActiveR() {return Active_diam/2.0/cm; }; 
    inline G4double GetActiveL() {return FielCageGap/cm; }; 
    inline G4double GetGasPressure(){return gas_pressure_;};
    inline G4double GetTemperature(){return temperature;};
    inline G4Region * GetRegion (){return rg;};
    inline G4LogicalVolume * GasLogic(){return gas_logic;}
    inline void SetMotherLab(G4LogicalVolume *mt){Mother=mt;};
    inline void SetOffset(G4double of){Offset=of;};


  
 private:
    DetectorMessenger* detectorMessenger;
    GasModelParameters* fGasModelParameters;
    G4bool checkOverlaps; // Check overlaps in the detector geometry if true
    G4double gas_pressure_; // pressure in the gas
    G4double temperature; // temperature of the gas

    G4double Lab_size;
    G4double chamber_diam   ;
    G4double chamber_length ;
    G4double chamber_thickn ;
    G4double SourceEn_offset ;
    G4double SourceEn_diam   ;
    G4double SourceEn_length ;
    G4double SourceEn_thickn ;
    G4double SourceEn_holedia ;
    G4ThreeVector vtx_;
    G4double Active_diam;
    G4double Active_length;
    G4double sc_yield_;
    G4double e_lifetime_;
    G4double ElGap_;
    G4double PMT1_Pos_;
    G4double PMT3_Pos_;
    G4ThreeVector vertex;
    pmt::PmtR7378A *pmt1_;
    pmt::PmtR7378A *pmt2_;
    G4LogicalVolume *Mother;

    G4double MgF2_window_thickness_;
    G4double MgF2_window_diam_;
    G4double pmt_hole_length_ ;
    G4double wndw_ring_stand_out_;
    G4double pedot_coating_thickness_;
    G4double optical_pad_thickness_;
    G4double pmt_base_diam_;
    G4double pmt_base_thickness_;
    G4bool  HideCollimator_;

    G4bool HideSourceHolder_;
    G4double max_step_size_;
    G4double ELyield_;

    G4double FielCageGap;
    G4double Offset;
    G4Region* rg;
    G4LogicalVolume* gas_logic;


};
#endif
