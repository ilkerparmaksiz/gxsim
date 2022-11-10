#ifndef GasModelParameters_hh
#define GasModelParameters_hh

#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include <map>

class DegradModel;
class GasModelParametersMessenger;
class DetectorConstruction;
class G4String;

class GasModelParameters{
	public:
	
	GasModelParameters();
	~GasModelParameters();
    
    /*Getters and Setters*/
    inline void SetThermalEnergy(G4double d){thermalE=d;}
    inline G4double GetThermalEnergy(){return thermalE;};

    void SetGasFile(G4String fn){gasFile = fn;}
    G4String GetGasFile(){return gasFile;};
    void SetFracXe(G4double fx){fracXe = fx;};
    G4double GetFracXe(){return fracXe;};
	
    private:
    GasModelParametersMessenger* fMessenger;
    G4double thermalE;
    G4String gasFile;
    G4double fracXe;
};

#endif
