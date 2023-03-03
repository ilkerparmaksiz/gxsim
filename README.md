[Eric's version of the following, with new application(s), like CRAB ....]

Geant4GarfieldDegradInterface
 ========================
 
 What is it?
 -----------
 A general example on how to interface Geant4 with Garfield++ and Degrad
 
 Documentation
 -------------
 Please check https://geant4garfielddegradinterface.readthedocs.io for a complete description of the implementation.
 
 Contributors
 ------------
 * Lennert De Keukeleere (lennert.dekeukeleere@kuleuven.be)
 * Dorothea Pfeiffer (Dorothea.Pfeiffer@cern.ch)


## INSTALLATION

In order to run the crab0 code, we need to install several software packages. The following instructions will help with that
```
NEST (forked from v2.3.12, 1 commit added for print statements)
Degrad (v3.15)
Garfield (master)
Geant4 v11.1.0 or greater
ROOT v6.24.04 works, but any recent stable version should be fine -- garfield needs it... sigh
GSL (just make sure it is setup and cmake should find it, I think...)
gxsim (main crab code)

```

I think for this to work you need a relatively up-to-date cmake v3 and gcc version. Let me know if you need the exact versioning if there are issues. 

---




Geant4 installation instructions

Download the geant4 tar file to the machine, you can scp it over to the packages area you just created
or use the following wget command to download it straight there:
```
wget https://geant4-data.web.cern.ch/geant4-data/releases/geant4-v11.1.0.tar.gz
```
when the copy has finished untar the geant4 tar file:
```
tar -xvf geant4-v11.0.0.tar.gz
cd geant4-v11.0.0
```

Make a geant4 build  and install directory and cd into the build:
```
mkdir build install
cd build
```

Ok with the prerequisits ready we can now install
```
cmake -DCMAKE_INSTALL_PREFIX=../install ../ -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_QT=OFF -DGEANT4_USE_GDML=ON -DCMAKE_CXX_FLAGS="-std=c++17"

make -j4
make install
```

---

ROOT installations

Sorry, gonna have to refer you this page, we all know how annoying ROOT is to build on any OS
https://root.cern/install/build_from_source/

Or if the machines have root compiled with c++17 instructions, just use that and life is easy.

---


Garfield Installation Instructions

k, now we are ready to build garfiled

Set some environmental variables to genat4 and root 
I am going to assume you now have all these set for the rest of the other packages
```
export G4INSTALL=/path/to/install;
export PATH=$G4INSTALL/bin:$PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;

cd $G4INSTALL/bin; source geant4.sh; cd -;

# Point to your root installation
export ROOTSYS=/path/to/root/6.24.04;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;
```
Now lets clone and install garfield
```
git clone https://gitlab.cern.ch/garfield/garfieldpp.git
cd garfieldpp
mkdir build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DWITH_EXAMPLES=OFF ../
make -j4
make install

# export the path to garfield
export GARFIELD_INSTALL=/Users/mistryk2/Packages/garfieldpp/install
export GARFIELD_HOME=/Users/mistryk2/Packages/garfieldpp/
export CMAKE_PREFIX_PATH=/Users/mistryk2/Packages/garfieldpp/install:$CMAKE_PREFIX_PATH
export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib
```

---

NEST Installation Instructions

Ok now we should try to install NEST
```
git clone https://github.com/kvjmistry/nest.git
git checkout crab0

cd nest
mkdir build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ../ -DG4=ON -DBUILD_ROOT=ON
make -j8
make install

#Set these environmental variables
export PATH=/path/to/nest/build/:$PATH;
export NEST_INCLUDE_DIRS=/path/to/nest/install/include/NEST

```

---

Degrad

The package can be downloaded from this website here:
https://degrad.web.cern.ch/degrad/
untar it in then compile with
```
gfortran -std=legacy -o Degrad degrad-3.15.f

# Set this path to your degrad folder
export DEGRAD_HOME=/Users/mistryk2/Packages/Degrad
```

---

The main repository can be downloaded from this github page and checking out the crab0 branch:
```
git clone https://github.com/kvjmistry/gxsim.git
git checkout crab0

# Set the path to the CRAB directory e.g.
export CRABPATH=/Users/mistryk2/Packages/GXeTPCSim/gxsim/CRAB/
```

```
# Make sure all the paths to the libraries are available
# here is an example setup script I have that does it all

# Set the path to the Geant4 Installation
export G4INSTALL=/Users/mistryk2/Packages/geant4-v11/geant4-v11.1.0/install;
export PATH=$G4INSTALL/bin:$PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;

cd $G4INSTALL/bin; source geant4.sh; cd -;


# Path to ROOT
export ROOTSYS=$(brew --cellar root)/6.24.04;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;

# Garfield
export GARFIELD_INSTALL=/Users/mistryk2/Packages/garfieldpp/install
export GARFIELD_HOME=/Users/mistryk2/Packages/garfieldpp/
export CMAKE_PREFIX_PATH=/Users/mistryk2/Packages/garfieldpp/install:$CMAKE_PREFIX_PATH
export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib

# DEGRAD
export DEGRAD_HOME=/Users/mistryk2/Packages/Degrad

# NEST
export NEST_INCLUDE_DIRS=/Users/mistryk2/Packages/NEST/install/include/NEST
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/mistryk2/Packages/NEST/install/lib

export CRABPATH=/Users/mistryk2/Packages/GXeTPCSim/gxsim/CRAB/

```

```
# Lets try and build it
mkdir build
cd build

# We need to point it to the gcem cmake file in nest for the cmake to be happy
cmake ../ -DCMAKE_PREFIX_PATH=/Users/mistryk2/Packages/NEST/install/lib/cmake/gcem
make -j4

# Add the crab exe to the path
export PATH=$CRABPATH/build:$PATH;
```

I have some data files that you need so you can run the comsol modes and particle generation files. You can copy them from this folder to your own CRAB/data folder on the BEBOP machine
```
/lcrc/project/NEXT/kmistry/software/gxsim/CRAB/data
```

---


Running CRAB

To run crab we can run the default macro in the macros folder, the seed number can just be any number greater than one.
```
# This runs an alpha at 5.3 MeV and a beta with energy spectrum up to 1.16 MeV
/path/to/build/CRAB macros/Alpha_e.mac <seed number>
```

Configuring the comsol modes in the macro file:
In the macro file you will see the following lines. Set the COMSOL path to the data folder where you copy over the comsol files. 
```
# Set this path to the data folder
/gasModelParameters/geometry/COMSOL_Path /Users/mistryk2/OneDrive - University of Texas at Arlington/Projects/CRAB/COMSOL/

# This uses the Garfield EL light model for generating S2. It reads e- trajectories from a file I made
/gasModelParameters/geometry/useEL_File false 

# This flag turns on comsol for the E field. 
/gasModelParameters/geometry/useComsol false # 
```

To control the event number you can set this number to what you want
```
/Action/SteppingAction/event_shift 0
```

