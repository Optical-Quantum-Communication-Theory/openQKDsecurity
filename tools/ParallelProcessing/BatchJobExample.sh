#!/usr/bin/env bash


# For this example, navigate to the openQKDSecurity main directory on command line then run
# ./tools/ParallelProcessing/BatchJobExample.sh

# If your using this on any sort of bash emulator on Windows, the paths automatically set bellow will not work with Matlab.
# Instead, manually write out the paths using the Windows path conventions instead. Don't forget escape character on \, \\.

# folder to start matlab / the working directory. This Default is where the script is run from.
WORKING_DIR=$(pwd);

# absolute or relative (from WORKING_DIR) to add aditional resources to Matlab's path. This Default is where the script is run from.
RESOURCE_DIR=$(pwd);

 # Name of the matlab preset/code to construct the QKDSolverInput object.
PRESET_NAME="BasicBB84WCPDecoyPreset";

# Absolute or relative (from WORKING_DIR) path to save output to.
FILE_OUT_DIR=$WORKING_DIR;

# Name you want individual job files (FILE_NAME_#.mat), log files (FILE_NAME_#.log), and the final compiled file(FILE_NAME.mat) to use.
FILE_NAME="protocolResults"; 

# How many jobs you want the preset to be split into. Can't exceed the total number of iterations for the QKDSolverInput.
# We recommend you keep this at most 1 less than the total number of physical cores on your system.
TOTAL_JOBS=5;

echo "Each MATLAB instance will be started in:" $WORKING_DIR
echo "Each MATLAB instance will load resources from:" $RESOURCE_DIR
echo "Each MATLAB instance will save files to:" $FILE_OUT_DIR

# Matlab is 1s based indexing, I'm sorry.
for ((JOB_IND = 1; JOB_IND <= TOTAL_JOBS; JOB_IND++)); do 

     echo "starting Job: ${JOB_IND}"
	 
	 # You may have to remove '-nojvm' if you encounter esoteric error messages.
     matlab -noFigureWindows -singleCompThread -nojvm -sd "${WORKING_DIR}" -batch "addpath(genpath('${RESOURCE_DIR}'));\
	 JobConstructor.runJobInLoop(${PRESET_NAME}, '${FILE_OUT_DIR}', '${FILE_NAME}', ${JOB_IND}, ${TOTAL_JOBS});" > \
	 "${FILE_NAME}_${JOB_IND}.log" &
done

wait

echo "Jobs Finished! Compiling Final Output..."

matlab -noFigureWindows -singleCompThread -nojvm -sd "${WORKING_DIR}" -batch "addpath(genpath('${RESOURCE_DIR}'));\
JobConstructor.compileJobs('${FILE_OUT_DIR}','${FILE_NAME}',true);"

echo "Finished Compiling"
