# Set the path to the Geant4 Installation
export G4INSTALL=/Users/mistryk2/Packages/geant4-v11/geant4-v11.1.0/install;
export PATH=$G4INSTALL/bin:$PATH;
export DYLD_LIBRARY_PATH=$G4INSTALL/lib:$DYLD_LIBRARY_PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;

cd $G4INSTALL/bin; source geant4.sh; cd -;


# Path to ROOT
export ROOTSYS=$(brew --cellar root)/6.24.04;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib;

# Garfield
export GARFIELD_INSTALL=/Users/mistryk2/Packages/garfieldpp/install
export GARFIELD_HOME=/Users/mistryk2/Packages/garfieldpp/
export CMAKE_PREFIX_PATH=/Users/mistryk2/Packages/garfieldpp/install:$CMAKE_PREFIX_PATH
export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$GARFIELD_INSTALL/lib



# DEGRAD
export DEGRAD_HOME=/Users/mistryk2/Packages/Degrad
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/homebrew/Cellar/gcc/11.3.0_2/lib;
export PATH=/opt/homebrew/bin:$PATH

# NEST
export NEST_INCLUDE_DIRS=/Users/mistryk2/Packages/NEST/install/include/NEST
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/mistryk2/Packages/NEST/install/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/mistryk2/Packages/NEST/install/lib;

export CRABPATH=/Users/mistryk2/Packages/GXeTPCSim/gxsim/CRAB/

# Add the crab exe to the path
export PATH=$CRABPATH/build:$PATH;
