//
// Created by argon on 2/15/24.
//

#ifndef CRAB_SIMPLEBOXGEO_HH
#define CRAB_SIMPLEBOXGEO_HH
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class SimpleBoxGeo {
public:
    SimpleBoxGeo();
    ~SimpleBoxGeo();
    void Construct();

    inline void CheckOverlaps(G4bool co){checkOverlaps=co;};
    inline void SetGasPressure(G4double d){gas_pressure_=d;};
    inline void SetTemperature(G4double d){temperature=d;};
    inline void SetGasLogic(G4LogicalVolume * gas){gas_logic=gas;};
    inline void SetYield(G4double  Y){Yield=Y;};
    //inline void SetMotherPhysical(G4VPhysicalVolume  *Ph){MotherPhysical=Ph;};
    inline G4double GetChamberR(){return chamber_diam/2.0/cm;};
    inline G4double GetChamberL(){return chamber_length/cm; };
    inline G4double GetActiveR() {return Active_diam/2.0/cm; };
    inline G4double GetActiveL() {return FielCageGap/cm; };
    inline G4double GetYield() {return Yield; };
    //inline G4VPhysicalVolume* GetMotherPhysical() {return MotherPhysical; };



    inline G4double GetGasPressure(){return gas_pressure_;};
    inline G4double GetTemperature(){return temperature;};
    inline G4LogicalVolume * GasLogic(){return gas_logic;}
    inline void SetMotherLab(G4LogicalVolume *mt){Mother=mt;};
private:
    G4LogicalVolume *Detector;

    G4bool checkOverlaps; // Check overlaps in the detector geometry if true
    G4double gas_pressure_; // pressure in the gas
    G4double temperature; // temperature of the gas

    G4double Lab_size;


    G4LogicalVolume *Mother;


    G4double chamber_diam   ;
    G4double chamber_length ;
    G4double Yield ;

    G4double Active_diam;
    G4double FielCageGap;

    G4double Offset;
    G4LogicalVolume* gas_logic;
    //G4VPhysicalVolume * MotherPhysical;
};


#endif //CRAB_SIMPLEBOXGEO_HH
