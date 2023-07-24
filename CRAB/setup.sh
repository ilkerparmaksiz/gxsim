# Set the path to the Geant4 Installation
export G4INSTALL=/home/argon/Programs/GEANT4/geant4-v11.1.1/install;
export PATH=$G4INSTALL/bin:$PATH;
export DYLD_LIBRARY_PATH=$G4INSTALL/lib:$DYLD_LIBRARY_PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;

cd $G4INSTALL/bin; source geant4.sh; cd -;


# Path to ROOT
export ROOTSYS=/home/argon/Programs/root_src/builddir;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib;

# Garfield
export GARFIELD_INSTALL=/home/argon/Programs/garfieldpp/install
export GARFIELD_HOME=/home/argon/Programs/garfieldpp/
export CMAKE_PREFIX_PATH=/home/argon/Programs/garfieldpp/install:$CMAKE_PREFIX_PATH
export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$GARFIELD_INSTALL/lib



# DEGRAD
export DEGRAD_HOME=/home/argon/Programs/Degrad
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/homebrew/Cellar/gcc/11.3.0_2/lib;
#export PATH=/opt/homebrew/bin:$PATH

# NEST
export NEST_INCLUDE_DIRS=/home/argon/Programs/nest/install/include/NEST
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/argon/Programs/nest/install/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/argon/Programs/nest/install/lib;

export gcem_DIR=/home/argon/Programs/nest/install/lib/cmake/gcem
export CRABPATH=/home/argon/Projects/Ilker/gxsim/CRAB

# Add the crab exe to the path
export PATH=$CRABPATH/build:$PATH;
