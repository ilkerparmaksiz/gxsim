#setenv PATH {$PATH}:/Applications/Emacs.app/Contents/MacOS:/Applications/CMake.app/Contents/bin:/usr/local/gfortran/bin

# Geant4 Path, edit G4Install path to where the main geant4 code folder your downloaded
export G4INSTALL=/Users/mistryk2/Packages//geant4-v10/geant4.10.06.p03;
export PATH=$G4INSTALL/bin:$PATH;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$G4INSTALL/lib;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4INSTALL/lib;

cd $G4INSTALL/bin;
source geant4.sh;
cd -;


# Path to ROOT
export ROOTSYS=$(brew --cellar root)/6.24.04;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib;

# Garfield
export GARFIELD_INSTALL=/Users/mistryk2/Packages/garfieldpp
export GARFIELD_HOME=/Users/mistryk2/Packages/garfieldpp
export CMAKE_PREFIX_PATH=/Users/mistryk2/Packages/garfieldpp:$CMAKE_PREFIX_PATH
export HEED_DATABASE=/Users/mistryk2/Packages/garfieldpp/Heed/heed++/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib

# DEGRAD
export DEGRAD_HOME=/Users/mistryk2/Packages/Degrad
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/homebrew/Cellar/gcc/11.3.0_2/lib;
export PATH=/opt/homebrew/bin:$PATH

# NEST
export NEST_INCLUDE_DIRS=/Users/mistryk2/Packages/NEST/install_g4v10/include/NEST
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/mistryk2/Packages/NEST/install_g4v10/lib


