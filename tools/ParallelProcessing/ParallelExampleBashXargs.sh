!/usr/bin/env bash

# MATLAB starting directory. Use START_DIR=$(pwd) to select the folder you
# called this script from.
START_DIR="path/to/project"

# absolute or relative (from START_DIR) path to add additional resources to
# MATLAB's path.
RESOURCE_DIR="path/to/resources"

# Name of the matlab preset/code to construct the QKDSolverInput object.
PRESET_NAME="BasicBB84WCPDecoyPreset"

# Absolute or relative (from START_DIR) path to save output files to.
FILE_OUT_DIR="path/to/save"

# Name you want individual job files (FILE_NAME_#.mat), log files
# (FILE_NAME_#.log), and the final compiled file(FILE_NAME.mat) to use.
FILE_NAME="BasicBB84WCPDecoyResults"

# Concatenate START_DIR with FILE_OUT_DIR for saving log files if
# FILE_OUT_DIR is a relative path.
if [[ "$FILE_OUT_DIR" == /* ]]; then
    LOG_DIR=$FILE_OUT_DIR # abs
else
    LOG_DIR="${START_DIR}/${FILE_OUT_DIR}" # rel
fi

# How many jobs you want the preset to be split into. Can't exceed the total
# number of iterations for the QKDSolverInput. We recommend a max of total
# physical cores -1 parallel jobs. For example, 6 cores -> 5 parallel jobs.
TOTAL_JOBS=7
TOTAL_PARALLEL_JOBS=5

echo
echo "Parallel Processing:" $PRESET_NAME
echo
echo "Start Directory:    " $START_DIR
echo "Resource Directory: " $RESOURCE_DIR
echo "Save Directory:     " $FILE_OUT_DIR
echo

# note every variable statement except {JOB_IND} is replaced before xargs
# replaces {JOB_IND}.
seq 1 "$TOTAL_JOBS" | xargs -P "$TOTAL_PARALLEL_JOBS" -I {JOB_IND} bash -c "
    echo \"starting Job: {JOB_IND}\"
    
    # You may have to remove '-nojvm' if you encounter esoteric error messages.
    matlab -noFigureWindows -singleCompThread -nojvm -sd \"${START_DIR}\" \
        -batch \"addpath(genpath('${RESOURCE_DIR}')); \
        JobConstructor.runJobInLoop(${PRESET_NAME}(), '${FILE_OUT_DIR}', \
        '${FILE_NAME}', {JOB_IND}, ${TOTAL_JOBS});\" > \
        \"${LOG_DIR}/${FILE_NAME}_{JOB_IND}.log\" 2>&1
"

echo "Jobs Finished! Compiling Final Output..."

# You may have to remove '-nojvm' if you encounter esoteric error messages.
matlab -noFigureWindows -singleCompThread -nojvm -sd "${START_DIR}" \
    -batch "addpath(genpath('${RESOURCE_DIR}')); \
    JobConstructor.compileJobs('${FILE_OUT_DIR}','${FILE_NAME}',true);"

echo "Finished Compiling"
