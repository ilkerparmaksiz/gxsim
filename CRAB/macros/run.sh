#!/bin/bash
## folder Name
folder=debug
## Set CRABPATH this is needed for the simulation
export CRABPATH=../${pwd}
echo "CRABPATH is $CRABPATH"
## if the folder does nt exist , this will create it
if [ ! -d "../$folder" ]
  then
  cd .. && mkdir "$folder" && cd "$folder" && cmake .. && make -j8
  cd ../macros
else
  cd .. && cd "$folder" && make -j8
  cd ../macros
fi
## if user passes int then interactive mode will be opened
## Run in interactive or batch mode
if [ "$1" = "int" ]; then
  echo "argument is $1"
  echo "Running Interactive Mode"
  ../"$folder"/CRAB
else
  echo "argument is $1"
  echo "Running Batch Mode"
  ../"$folder"/CRAB Single_OpticalPhoton.mac 100
fi

