# MVP power shell batch job example 

# MATLAB starting directory. Use $START_DIR=((Get-Location).Path) to select
# the folder you called this script from.
$START_DIR="path/to/project"

# absolute or relative (from START_DIR) path to add additional resources to
# MATLAB's path.
$RESOURCE_DIR="path/to/resources"

# Name of the matlab preset/code to construct the QKDSolverInput object.
$PRESET_NAME="BasicBB84WCPDecoyPreset"

# Absolute or relative (from START_DIR) path to save output files to.
$FILE_OUT_DIR="path/to/save"

# Name you want individual job files (FILE_NAME_#.mat), log files
# (FILE_NAME_#.log), and the final compiled file(FILE_NAME.mat) to use.
$FILE_NAME="BasicBB84WCPDecoyResults"

if ([System.IO.Path]::IsPathRooted($FILE_OUT_DIR)){
    # abs path
    $LOG_DIR="$FILE_OUT_DIR"
} else {
    # rel path
    $LOG_DIR=(Join-Path -Path "$START_DIR" -ChildPath "$FILE_OUT_DIR")
}

# How many jobs you want the preset to be split into. Can't exceed the total
# number of iterations for the QKDSolverInput. We recommend a max of total
# physical cores -1 parallel jobs. For example, 6 cores -> 5 parallel jobs.
$TOTAL_JOBS=7
$TOTAL_PARALLEL_JOBS=5

Write-Output ""
Write-Output "Parallel Processing: $PRESET_NAME"
Write-Output ""
Write-Output "Start Directory:     $START_DIR"
Write-Output "Resource Directory:  $RESOURCE_DIR"
Write-Output "Save Directory:      $FILE_OUT_DIR"
Write-Output ""

# Matlab is 1s based indexing, I'm sorry.
(1..$TOTAL_JOBS) | ForEach-Object -Parallel {

    Write-Output "Starting Job: $_"

    $MATLAB_INPUT =("addpath(genpath('{0}')); `
        JobConstructor.runJobInLoop({1}(), '{2}', '{3}', {4}, {5});" `
        -f $using:RESOURCE_DIR, $using:PRESET_NAME, `
        $using:FILE_OUT_DIR, $using:FILE_NAME, $_, $using:TOTAL_JOBS)

    $LOG_PATH_AND_NAME = (Join-Path -Path "$using:LOG_DIR" -ChildPath `
        (("{0}_{1}.log" -f $using:FILE_NAME, $_)))

    # You may have to remove '-nojvm' if you encounter esoteric error messages.
    matlab -noFigureWindows -singleCompThread -nojvm -sd $using:START_DIR `
        -batch $MATLAB_INPUT `
        > "$LOG_PATH_AND_NAME" 2>&1

} -ThrottleLimit $TOTAL_PARALLEL_JOBS

Write-Output "Jobs Finished! Compiling Final Output..."

# You may have to remove '-nojvm' if you encounter esoteric error messages.
matlab -noFigureWindows -singleCompThread -nojvm -sd $START_DIR `
    -batch ("addpath(genpath('{0}')); JobConstructor.compileJobs('{1}','{2}', true);" `
    -f $RESOURCE_DIR, $FILE_OUT_DIR, $FILE_NAME)

Write-Output "Finished Compiling"