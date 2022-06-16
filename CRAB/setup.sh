alias emacs Emacs

setenv PATH {$PATH}:/Applications/Emacs.app/Contents/MacOS:/Applications/CMake.app/Contents/bin:/usr/local/gfortran/bin

cd /Users/chur558/geant4.10.07.02-install/bin
source /Users/chur558/geant4.10.07.02-install/bin/geant4.csh
cd -
source /Users/chur558/root_build/bin/thisroot.csh

setenv MARLEY /Users/chur558/marley
setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/Users/chur558/marley/build:/usr/local/gfortran/lib
setenv DYLD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/Users/chur558/marley/build:/usr/local/gfortran/lib

setenv GARFIELD_HOME /Users/chur558/we22903/garfieldpp
setenv HEED_DATABASE /Users/chur558/we22903/garfieldpp/Heed/heed++/database
setenv DEGRAD_HOME /Users/chur558/degrad

setenv NEST_INCLUDE_DIRS /Users/chur558/nest//install/include/NEST

source ~/venvs/bin/activate.csh

# to get proper tab completion in csh, it seems.
set filec
set autolist
