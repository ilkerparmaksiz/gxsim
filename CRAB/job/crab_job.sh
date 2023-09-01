#!/bin/bash
#SBATCH -J CRAB # A single job name for the array
#SBATCH -c 1 # Number of cores
#SBATCH -p node2
#SBATCH --mem 30000 # Memory request (6Gb)
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o CRAB_%A_%a.out # Standard output
#SBATCH -e CRAB_%A_%a.err # Standard error

start=`date +%s`

# Set the configurable variables
JOBNAME="Alpha"
TYPE="CRAB"
N_EVENTS=5

# Create the directory
source "/home/argon/Projects/Ilker/gxsim/CRAB/macros/run.sh test"
cd /media/argon/Data/CRAB/Sim
mkdir -p $JOBNAME/$TYPE/jobid_"${SLURM_ARRAY_TASK_ID}"
cd $JOBNAME/$TYPE/jobid_"${SLURM_ARRAY_TASK_ID}"

# Copy the macro file
cp $CRABPATH/macros/Single_alpha.mac .

# Setup nexus and run
echo "Setting Up Code" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
#source /home/argon/Projects/Krishan/gxsim/CRAB/setup_cluster.sh

# Calculate the unique seed number	
SEED=$((${N_EVENTS}*(${SLURM_ARRAY_TASK_ID} - 1) + ${N_EVENTS}))
echo "The seed number is: ${SEED}" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt

# Replace the number of events in the file as well as the event index
sed -i "s#.*event_shift.*#/Action/SteppingAction/event_shift ${SEED}#" run1.mac
sed -i "s#.*beamOn.*#/run/beamOn ${N_EVENTS}#" run1.mac

# NEXUS
echo "Running GXeSim" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
$CRABPATH/build/CRAB Single_alpha.mac ${SEED} 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt

echo; echo; echo;

echo "FINISHED....EXITING" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
