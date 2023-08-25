// ----------------------------------------------------------------------------
// nexus | OpticalMaterialProperties.cc
//
// Optical properties of relevant materials.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticalMaterialProperties.hh"
#include "XenonProperties.hh"
#include "SellmeierEquation.hh"

#include <G4MaterialPropertiesTable.hh>

#include <assert.h>

using namespace CLHEP;

namespace opticalprops {
  /// Vacuum ///
  G4MaterialPropertiesTable* Vacuum()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> photEnergy = {optPhotMinE_, optPhotMaxE_};

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex = {1., 1.};
    mpt->AddProperty("RINDEX", photEnergy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> absLength = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", photEnergy, absLength);

    return mpt;
  }

  /// Fused Silica ///
  G4MaterialPropertiesTable* FusedSilica()
  {
    // Optical properties of Suprasil 311/312(c) synthetic fused silica.
    // Obtained from http://heraeus-quarzglas.com

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
    // for fused silica is valid only in that range
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    // The following values for the refractive index have been calculated
    // using Sellmeier's equation:
    //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
    //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
    //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
    // with wavelength \lambda in micrometers and
    //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
    //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

    G4double B_1 = 4.73e-1;
    G4double B_2 = 6.31e-1;
    G4double B_3 = 9.06e-1;
    G4double C_1 = 1.30e-2;
    G4double C_2 = 4.13e-3;
    G4double C_3 = 9.88e+1;

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
        + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
        + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
      rIndex.push_back(sqrt(n2));
      // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,  6.46499 * eV,
      6.54000 * eV,  6.59490 * eV,  6.64000 * eV,  6.72714 * eV,
      6.73828 * eV,  6.75000 * eV,  6.82104 * eV,  6.86000 * eV,
      6.88000 * eV,  6.89000 * eV,  7.00000 * eV,  7.01000 * eV,
      7.01797 * eV,  7.05000 * eV,  7.08000 * eV,  7.08482 * eV,
      7.30000 * eV,  7.36000 * eV,  7.40000 * eV,  7.48000 * eV,
      7.52000 * eV,  7.58000 * eV,  7.67440 * eV,  7.76000 * eV,
      7.89000 * eV,  7.93000 * eV,  8.00000 * eV,
      optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      200.0 * cm,   200.0 * cm,  90.0 * cm,  45.0 * cm,
      45.0 * cm,    30.0 * cm,  24.0 * cm,  21.0 * cm,
      20.0 * cm,    19.0 * cm,  16.0 * cm,  14.0 * cm,
      13.0 * cm,     8.5 * cm,   8.0 * cm,   6.0 * cm,
       1.5 * cm,     1.2 * cm,   1.0 * cm,   .65 * cm,
        .4 * cm,     .37 * cm,   .32 * cm,   .28 * cm,
        .22 * cm,    .215 * cm,  .00005*cm,
      .00005* cm
    };

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }


  // Absorbtion Length https://docs.google.com/spreadsheets/d/1MtjtvoWGBKZvObom2R21wITAbv8SQGjoolY9mlpAGcs/edit?usp=sharing
  G4MaterialPropertiesTable * MgF2(){
    
    G4MaterialPropertiesTable *mpt =new G4MaterialPropertiesTable();
        std::vector<G4double> RIndex;
        // REFRACTIVE INDEX
        //https://refractiveindex.info/?shelf=main&book=MgF2&page=Dodge-o
        G4double um2 = micrometer*micrometer;
        G4double B[3] = {0.48755108, 0.39875031	, 2.3120353};
        G4double C[3] = {0.001882178 * um2, 0.008951888 * um2, 566.13559 * um2};
        SellmeierEquation seq(B, C);

        const G4int ri_entries = 100;
        G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

        std::vector<G4double> ri_energy;
        for (int i=0; i<ri_entries; i++)
        {
              ri_energy.push_back(optPhotMinE_ + i * eWidth);
              RIndex.push_back(seq.RefractiveIndex(h_Planck*c_light/ri_energy[i]));
            //G4cout << "* MgF2 rIndex:  " << std::setw(5)
                 // << (h_Planck*c_light/ri_energy[i])/nm << " nm -> " << RIndex[i] << G4endl;
        }
          mpt->AddProperty("RINDEX", ri_energy, RIndex);
        // AbsLength
        std::vector<G4double> AbsEnergy;
        std::vector<G4double> AbsLength;
    
    return mpt;
}


  /// Sapphire ///
  G4MaterialPropertiesTable* Sapphire()
  {
    // Input data: Sellmeier equation coeficients extracted from:
    // https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire
    //C[i] coeficients at line 362 are squared.

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {1.4313493, 0.65054713, 5.3414021};
    G4double C[3] = {0.0052799261 * um2, 0.0142382647 * um2, 325.017834 * um2};
    SellmeierEquation seq(B, C);

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(seq.RefractiveIndex(h_Planck*c_light/ri_energy[i]));
      //G4cout << "* Sapphire rIndex:  " << std::setw(5)
      //       << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_, 0.900 * eV,
      1.000 * eV,   1.296 * eV,  1.683 * eV,  2.075 * eV,
      2.585 * eV,   3.088 * eV,  3.709 * eV,  4.385 * eV,
      4.972 * eV,   5.608 * eV,  6.066 * eV,  6.426 * eV,
      6.806 * eV,   7.135 * eV,  7.401 * eV,  7.637 * eV,
      7.880 * eV,   8.217 * eV
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      3455.0  * mm,  3455.0  * mm,  3455.0  * mm,  3455.0  * mm,
      3455.0  * mm,  3140.98 * mm,  2283.30 * mm,  1742.11 * mm,
      437.06 * mm,   219.24 * mm,  117.773 * mm,   80.560 * mm,
      48.071 * mm,   28.805 * mm,   17.880 * mm,   11.567 * mm,
        7.718 * mm,    4.995 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



   /// Gaseous Xenon ///
  G4MaterialPropertiesTable* GXe(G4double pressure,
                                G4double temperature,
                                G4int    sc_yield,
                                G4double e_lifetime)
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    G4double density = GXeDensity(pressure);
    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(XenonRefractiveIndex(ri_energy[i], density));
      // G4cout << "* GXe rIndex:  " << std::setw(7)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex, ri_entries);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // EMISSION SPECTRUM
    // Sampling from ~150 nm to 200 nm <----> from 6.20625 eV to 8.20625 eV
    const G4int sc_entries = 200;
    std::vector<G4double> sc_energy;
    for (int i=0; i<sc_entries; i++){
      sc_energy.push_back(6.20625 * eV + 0.01 * i * eV);
    }
    std::vector<G4double> intensity;
    for (G4int i=0; i<sc_entries; i++) {
      intensity.push_back(GXeScintillation(sc_energy[i], pressure));
    }
    //for (int i=0; i<sc_entries; i++) {
    //  G4cout << "* GXe Scint:  " << std::setw(7) << sc_energy[i]/eV
    //         << " eV -> " << intensity[i] << G4endl;
    //}
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);
    mpt->AddProperty("ELSPECTRUM"             , sc_energy, intensity, 1);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", sc_yield);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   4.5  * ns);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   100. * ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .1);
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .9);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime, 1);

    return mpt;
  }
  
  // Stainles Steel Optical Properties Table
  G4MaterialPropertiesTable * STEEL()
  {
      G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

      // REFLECTIVITY
      std::vector<G4double> ENERGIES = {
              optPhotMinE_,7.20*eV, 7.29 * eV,  optPhotMaxE_
      };
      //std::vector<G4double> REFLECTIVITY = { 0.20,0.20,0.20,0.20};
      std::vector<G4double> REFLECTIVITY = { 0.2,0.2,0.2,0.2};
      // std::vector<G4double> REFLECTIVITY = { 0.00,0.00,0.00};

      // REFLEXION BEHAVIOR
      std::vector<G4double> ENERGIES_2    = {optPhotMinE_, optPhotMaxE_};
      // Specular reflection about the normal to a microfacet.
      // Such a vector is chosen according to a gaussian distribution with
      // sigma = SigmaAlhpa (in rad) and centered in the average normal.
      std::vector<G4double> specularlobe  = {0., 0.};
      // specular reflection about the average normal
      std::vector<G4double> specularspike = {0., 0.};
      // 180 degrees reflection.
      std::vector<G4double> backscatter   = {0., 0.};
      // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

      mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe);
      mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike);
      mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter);
      mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);
      return mpt;
  }
  G4MaterialPropertiesTable * PerfectDetector(){
      G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

      // REFLECTIVITY
      std::vector<G4double> ENERGIES = {
              optPhotMinE_,7.20*eV, 7.29 * eV,  optPhotMaxE_
      };
      std::vector<G4double> REFLECTIVITY = { 0,0,0,0};
      std::vector<G4double> EFFICIENCY = { 1,1,1,1};

      mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);
      mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY);
      return mpt;
  }



}
