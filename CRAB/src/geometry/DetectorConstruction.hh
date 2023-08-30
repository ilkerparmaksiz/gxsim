#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh 1

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
#include "SampleFromSurface.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UniformMagField;


using namespace std;
/*! \class  DetectorConstruction*/
/*! \brief class derived from G4VUserDetectorConstruction*/

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
    DetectorConstruction(GasModelParameters*);
    virtual ~DetectorConstruction();

    // Mandatory methods
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    virtual void AssignVisuals();

    //Setters for the dimensions and environment variables of the setup
    inline void CheckOverlaps(G4bool co){checkOverlaps=co;};
    inline void SetGasPressure(G4double d){gas_pressure_=d;};
    inline void SetTemperature(G4double d){temperature=d;};

    inline G4double GetChamberR(){return chamber_diam/2.0/cm;};
    inline G4double GetChamberL(){return chamber_length/cm; }; 
    inline G4double GetActiveR() {return Active_diam/2.0/cm; }; 
    inline G4double GetActiveL() {return FielCageGap/cm; }; 
    inline G4double GetGasPressure(){return gas_pressure_;};
    inline G4double GetTemperature(){return temperature;};
    inline G4double GetELPosition(){return fEL_Pos;};
    inline G4double GetOffset(){return fOffset/cm;};
    inline G4double SetELPosition(G4double k){return fEL_Pos=k;};

    // Garfield  and Handling Multi-Threading
    inline  std::shared_ptr<Garfield::ComponentComsol> getComsolComponent (){return fcomsol;};
    inline  Garfield::MediumMagboltz* getMediumMagBoltz(){return fMediumMagboltz;};
    inline  Garfield::AvalancheMC * getAvalanceMC (){return fAvalancheMC;};
    inline  Garfield::Sensor*  getSensor (){return fSensor;};
    inline  std::vector<std::vector<std::vector<G4double>>>  getELProfiles (){return EL_profiles;};
    inline  std::vector<G4double> getELEvents (){return EL_events;};
    //const util::SampleFromSurface * getSamples() {return Sampler;};
    void GarfieldInitialization();
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
    G4double Gas_diam   ;
    G4double Gas_length ;
    std::shared_ptr<Garfield::ComponentComsol> fcomsol ;
    G4ThreeVector vtx_;
    G4double Active_diam;
    G4double Active_length;
    G4double sc_yield_;
    G4double e_lifetime_;
    G4double ElGap_;
    G4double fEL_Pos;
    G4double fOffset;
    G4double counter;
    G4ThreeVector vertex;

    // Variable to store the EL timing profiles to sample from
    // <event> <photon> <x,y,z,t of photon>
    std::vector<std::vector<std::vector<G4double>>> EL_profiles;

    // Vector consisting of the event numbers from the simulated Garfield using
    // COMSOL geometry
    std::vector<G4double> EL_events;

    G4double MgF2_window_thickness_;
    G4double MgF2_window_diam_;


    G4bool HideSourceHolder_;
    G4double max_step_size_;
    G4double ELyield_;

    G4double FielCageGap;

    G4LogicalVolume* gas_logic;
    std::shared_ptr<util::SampleFromSurface>Sampler;

    //Garfield
    G4String gasFile;
    G4String ionMobFile;
    Garfield::MediumMagboltz* fMediumMagboltz;
    Garfield::AvalancheMC * fAvalancheMC;

    Garfield::Sensor* fSensor;
    filehandler::FileHandling FileHandler;

};
#endif
