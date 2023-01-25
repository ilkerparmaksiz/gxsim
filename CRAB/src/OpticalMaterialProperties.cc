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

    G4double photEnergy[] = {optPhotMinE_, optPhotMaxE_};
    G4double nEntries = sizeof(photEnergy) / sizeof(G4double);

    // REFRACTIVE INDEX
    G4double rIndex[] = {1., 1.};
    assert(sizeof(rIndex) == sizeof(photEnergy));
    mpt->AddProperty("RINDEX", photEnergy, rIndex, nEntries);

    // ABSORPTION LENGTH
    G4double absLength[] = {noAbsLength_, noAbsLength_};
    assert(sizeof(absLength) == sizeof(photEnergy));
    mpt->AddProperty("ABSLENGTH", photEnergy, absLength, nEntries);

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

    G4double ri_energy[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      ri_energy[i] = optPhotMinE_ + i * eWidth;
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

    G4double rIndex[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
        + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
        + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
      rIndex[i] = sqrt(n2);
      // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }

    G4double nEntries_ri = sizeof(ri_energy) / sizeof(G4double);
    mpt->AddProperty("RINDEX", ri_energy, rIndex, nEntries_ri);

    // ABSORPTION LENGTH
    G4double abs_energy[] = {
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

    G4double absLength[] = {
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

    G4double nEntries_abs = sizeof(abs_energy) / sizeof(G4double);
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength, nEntries_abs);

    return mpt;
  }


  // Absorbtion Length https://docs.google.com/spreadsheets/d/1MtjtvoWGBKZvObom2R21wITAbv8SQGjoolY9mlpAGcs/edit?usp=sharing
  G4MaterialPropertiesTable * MgF2(){
    
    G4MaterialPropertiesTable *mpt =new G4MaterialPropertiesTable();
    
    
    // REFRACTIVE INDEX
    //https://refractiveindex.info/?shelf=main&book=MgF2&page=Dodge-o
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {0.48755108, 0.39875031	, 2.3120353};
    G4double C[3] = {0.001882178 * um2, 0.008951888 * um2, 566.13559 * um2};
    SellmeierEquation seq(B, C);

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    G4double ri_energy[ri_entries];
    G4double RIndex[ri_entries];
    G4double MgF2Reflectivity[ri_entries];
    G4double MgF2Efficiency[ri_entries];
    for (int i=0; i<ri_entries; i++)
    {
          ri_energy[i] = optPhotMinE_ + i * eWidth;
          RIndex[i] = seq.RefractiveIndex(h_Planck*c_light/ri_energy[i]);
          MgF2Reflectivity[i] = 0.;
          MgF2Efficiency[i] = 1.0;
    }


    G4double nEntries = sizeof(ri_energy) / sizeof(G4double);

    mpt->AddProperty("RINDEX", ri_energy, RIndex, nEntries);
    mpt->AddProperty("REFLECTIVITY",ri_energy,MgF2Reflectivity,nEntries);
    mpt->AddProperty("EFFICIENCY",ri_energy,MgF2Efficiency,nEntries);
    
    return mpt;
}


  /// Sapphire ///
  G4MaterialPropertiesTable* Sapphire()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {1.4313493, 0.65054713, 5.3414021};
    G4double C[3] = {0.0052799261 * um2, 0.0142382647 * um2, 325.017834 * um2};
    SellmeierEquation seq(B, C);

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    G4double ri_energy[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      ri_energy[i] = optPhotMinE_ + i * eWidth;
    }

    G4double rIndex[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      rIndex[i] = seq.RefractiveIndex(h_Planck*c_light/ri_energy[i]);
      //G4cout << "* Sapphire rIndex:  " << std::setw(5)
      //       << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    assert(sizeof(rIndex) == sizeof(ri_energy));
    mpt->AddProperty("RINDEX", ri_energy, rIndex, ri_entries);

    // ABSORPTION LENGTH
    G4double abs_energy[] = {
      optPhotMinE_, 0.900 * eV,
      1.000 * eV,   1.296 * eV,  1.683 * eV,  2.075 * eV,
      2.585 * eV,   3.088 * eV,  3.709 * eV,  4.385 * eV,
      4.972 * eV,   5.608 * eV,  6.066 * eV,  6.426 * eV,
      6.806 * eV,   7.135 * eV,  7.401 * eV,  7.637 * eV,
      7.880 * eV,   8.217 * eV
    };
    const G4int abs_entries = sizeof(abs_energy) / sizeof(G4double);

    G4double absLength[] = {
      noAbsLength_, noAbsLength_,
      3455.0  * mm,  3455.0  * mm,  3455.0  * mm,  3455.0  * mm,
      3455.0  * mm,  3140.98 * mm,  2283.30 * mm,  1742.11 * mm,
      437.06 * mm,   219.24 * mm,  117.773 * mm,   80.560 * mm,
      48.071 * mm,   28.805 * mm,   17.880 * mm,   11.567 * mm,
        7.718 * mm,    4.995 * mm
    };
    assert(sizeof(absLength) == sizeof(abs_energy));
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength, abs_entries);

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

    G4double ri_energy[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      ri_energy[i] = optPhotMinE_ + i * eWidth;
    }

    G4double density = GXeDensity(pressure);
    G4double rIndex[ri_entries];
    for (int i=0; i<ri_entries; i++) {
      rIndex[i] = XenonRefractiveIndex(ri_energy[i], density);
      // G4cout << "* GXe rIndex:  " << std::setw(7)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    assert(sizeof(rIndex) == sizeof(ri_energy));
    mpt->AddProperty("RINDEX", ri_energy, rIndex, ri_entries);

    // ABSORPTION LENGTH
    G4double abs_energy[] = {optPhotMinE_, optPhotMaxE_};
    G4double absLength[]  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength, 2);

    // EMISSION SPECTRUM
    // Sampling from ~150 nm to 200 nm <----> from 6.20625 eV to 8.20625 eV
    const G4int sc_entries = 200;
    G4double sc_energy[sc_entries];
    for (int i=0; i<sc_entries; i++){
      sc_energy[i] = 6.20625 * eV + 0.01 * i * eV;
    }
    G4double intensity[sc_entries];
    // XenonScintillation(sc_entries, sc_energy, intensity, pressure);

    for (G4int i=0; i<sc_entries; i++) {
      intensity[i] = GXeScintillation(sc_energy[i], pressure);
    }

    //for (int i=0; i<sc_entries; i++) {
    //  G4cout << "* GXe Scint:  " << std::setw(7) << sc_energy[i]/eV
    //         << " eV -> " << intensity[i] << G4endl;
    //}
    mpt->AddProperty("FASTCOMPONENT", sc_energy, intensity, sc_entries);
    mpt->AddProperty("ELSPECTRUM",    sc_energy, intensity, sc_entries);
    mpt->AddProperty("SLOWCOMPONENT", sc_energy, intensity, sc_entries);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", sc_yield);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("FASTTIMECONSTANT",   4.5  * ns);
    mpt->AddConstProperty("SLOWTIMECONSTANT",   100. * ns);
    mpt->AddConstProperty("YIELDRATIO",         .1);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime);

    return mpt;
  }
  
  // Stainles Steel Optical Properties Table
  G4MaterialPropertiesTable * STEEL()
  {
      G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

      // REFLECTIVITY
      G4double ENERGIES[] = {
              optPhotMinE_, 7.29 * eV,  optPhotMaxE_
      };
      G4double REFLECTIVITY[] = {
              0.20,0.20,0.20
      };

      // REFLEXION BEHAVIOR
      G4double ENERGIES_2[]    = {optPhotMinE_, optPhotMaxE_};
      // Specular reflection about the normal to a microfacet.
      // Such a vector is chosen according to a gaussian distribution with
      // sigma = SigmaAlhpa (in rad) and centered in the average normal.
      G4double specularlobe[]  = {0., 0.};
      // specular reflection about the average normal
      G4double specularspike[] = {0., 0.};
      // 180 degrees reflection.
      G4double backscatter[]   = {0., 0.};
      // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

      G4double nEntries = sizeof(ENERGIES_2) / sizeof(G4double);

      mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe,nEntries);
      mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike,nEntries);
      mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter,nEntries);
      mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY,nEntries);
      return mpt;
  }


}
