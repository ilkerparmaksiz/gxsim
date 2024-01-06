#!/bin/bash
## folder Name
folder=build
num_photons=100
## Set CRABPATH this is needed for the simulation
export CRABPATH=${PWD%macros}
alias mc='(cd "${CRABPATH}macros")'
#echo "CRABPATH is $CRABPATH"
## if the folder does nt exist , this will create it
build(){
  if [ ! -d "../$folder" ]
    then
    cd .. && mkdir "$folder" && cd "$folder" && cmake .. && make -j8
    cd ../macros
  else
    cd .. && cd "$folder" && make -j8
    cd ../macros
  fi
}
#Run the build

## if user passes int then interactive mode will be opened
## Run in interactive or batch mode
if [ "$1" = "int" ]; then
  echo "argument is $1"
  echo "Running Interactive Mode"
  build
  ../"$folder"/CRAB
elif [ "$1" = "clean" ]; then
  echo "argument is $1"
  rm -rf ${CRABPATH}"$folder"
  build
  #echo "Running Interactive Mode"
  #../"$folder"/CRAB Single_OpticalPhoton.mac 100
elif [ "$1" = "mthread" ]; then
  echo "argument is $1"
  build
  echo "Running MultiThread Batch Mode"
   ../"$folder"/CRAB Single_OpticalPhoton.mac 100 30
elif [ "$1" = "malpha" ]; then
  echo "argument is $1"
  build
  echo "Running MultiThread Batch Mode"
   ../"$folder"/CRAB Single_alpha.mac 100 $2
elif [ "$1" = "alpha" ]; then
  echo "argument is $1"
  build
  echo "Running MultiThread Batch Mode"
   ../"$folder"/CRAB Single_alpha.mac 100
elif [ "$1" = "test" ]; then
  echo "CRABPATH is $CRABPATH";
  echo "Folder Path is $folder";  
elif [ "$1" = "optical" ]; then
  echo "argument is $1"
  echo "Running Batch Mode"
  build
  ../"$folder"/CRAB Single_OpticalPhoton.mac 100
else
  echo "argument is $1"
  echo "Running Batch Mode"
  build
  ../"$folder"/CRAB Single_OpticalPhoton.mac 100	
fi


