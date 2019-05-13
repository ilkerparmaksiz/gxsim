#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4SystemOfUnits.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*! \class DetectorMessenger*/
/*! class derived from G4UImessenger*/
/*! List of available commands*/
/*!/MWGPC/geometry/SetGasPressure*/
/*!/MWGPC/geometry/SetNrUpperPlanes*/
/*!/MWGPC/geometry/SetMaxStep*/
/*!/MWGPC/geometry/SetBField*/
/*!/MWGPC/geometry/EMField_version */
/*!/MWGPC/geometry/ConstructWires */
/*!/MWGPC/geometry/ConstructSlits5Vertical */
/*!/MWGPC/geometry/ConstructSlits3Vertical */
/*!/MWGPC/geometry/ConstructSlitVertical */
/*!/MWGPC/geometry/ConstructSlitHorizontal */
/*!/MWGPC/geometry/buildCells*/
/*!/MWGPC/geometry/BuildUpperScint*/
/*!/MWGPC/geometry/BuildLowerScint*/
/*!/MWGPC/geometry/update */

class DetectorMessenger : public G4UImessenger {
 public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DetectorConstruction* detector;

  G4UIdirectory* miniDir;      ///<\brief /MWGPC/
  G4UIdirectory* geometryDir;  ///<\brief /MWGPC/geometry/

  G4UIcmdWithADoubleAndUnit* setGasPressCmd;
  G4UIcmdWithAnInteger* setNumHexesCmd;
  G4UIcmdWithAString* setupNameCmd;
    
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
