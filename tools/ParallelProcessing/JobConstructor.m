classdef JobConstructor
    % JobConstructor tool to assist parallel processing outside of Matlab.
    % Because packages like CVX use global variables, we cannot reliably
    % use parfor or other parallel features from the parallel compute
    % toolbox. To get around this we can spawn multiple headless instances
    % of Matlab and get each to compute the key rates for a subset of the
    % points. This class is designed to help make that task a little easier
    % by providing functions to help break down a QKDSolverInput preset
    % into multiple slices, save individual results, and compile them into
    % a final output file.
    %
    % You can find an example in "BatchJobExample.sh" for Linux.
    %
    % See also: MainIteration, QKDSolverInput

    methods (Static)

        function runJob(qkdInput,basePath,baseFileName,jobIndex,scanIndexList)
            % Runs the QKDSolverInput for only the points specified in
            % scanIndexList. Then, it saves the results to the specified
            % folder under the  file name 'filename_jobIndex.mat'. The
            % scanIndexList is also included in the the file. Use
            % JobConstructor.compileJobs to merge all the runs together.
            %
            %
            % Inputs:
            % * qkdInput: A full QKDSolverInput with all required modules
            %   passed to MainIteration.
            % * basePath: String for the absolute or relative path to the
            %   folder where result files should be saved. If you construct
            %   this path in Matlab, we recommend using the function
            %   FULLFILE for OS agnostic file paths.
            % * baseFileName: The base file name used for output. The file
            %   will be named 'baseFileName_jobIndex.mat' to identify it.
            % * jobIndex: Positive 64 bit unsigned integer to distinguish
            %   this job from others.
            % * scanIndexList: Nonempty array of 64 bit unsigned integers
            %   that dictate which points from the scan parameters
            %   (compiled into linear indexing) should be run. The list
            %   must be a subset of 1:qkdInput.totalIterations.
            %
            % See also: MainIteration, QKDSolverInput, fullfile, JobConstructor.compileJobs
            arguments
                qkdInput (1,1) QKDSolverInput
                basePath (1,1) string {mustBeFolder}
                baseFileName (1,1) string {checkBaseFileName}
                jobIndex (1,1) uint64 {mustBePositive}
                scanIndexList (:,1) uint64 {mustBeNonempty,mustBePositive, ...
                    mustBeInTotalScanParamRange(scanIndexList,qkdInput)}
            end

            results = MainIteration(qkdInput,scanIndexList);

            fileName = sprintf("%s_%d.mat",baseFileName,jobIndex);
            fileName = fullfile(basePath,fileName);

            save(fileName,"results","qkdInput","scanIndexList");
        end



        %% running jobs
        function runJobStartEnd(qkdInput,basePath,baseFileName,jobIndex,startIndex,endIndex)
            % This is a small wrapper for JobConstructor.runJob which runs
            % a QKDSolverInput for all indexes in the range
            % startIndex:endIndex instead of passing in an array of
            % indexes.
            %
            % Inputs:
            % * qkdInput: A full QKDSolverInput with all required modules
            %   passed to MainIteration.
            % * basePath: string for the absolute or relative path to the
            %   folder where result files should be saved. If you construct
            %   this path in Matlab, we recommend using the function
            %   FULLFILE for OS agnostic file paths.
            % * baseFileName: The base file name used for output. The file
            %   will be named 'baseFileName_jobIndex.mat' to identify it.
            % * jobIndex: Positive 64 bit unsigned integer to distinguish
            %   this job from others.
            % * startIndex:
            % * endIndex:
            %
            % See also: MainIteration, QKDSolverInput, fullfile,
            % JobConstructor.runJob
            arguments
                qkdInput (1,1) QKDSolverInput
                basePath (1,1) string {mustBeFolder}
                baseFileName (1,1) string {checkBaseFileName}
                jobIndex (1,1) uint64 {mustBePositive}
                startIndex (1,1) uint64 {mustBePositive}
                endIndex (1,1) uint64 {mustBeGreaterThanOrEqual(endIndex,startIndex)}
            end
            JobConstructor.runJob(qkdInput,basePath,baseFileName,jobIndex,startIndex:endIndex);
        end



        function runJobInLoop(qkdInput,basePath,baseFileName,jobIndex,numJobs)
            % This is a wrapper for JobConstructor.runJob which equally
            % divides the QKDSolverInput into the specified number of jobs.
            % Call this function in a loop (for example, in a bash script)
            % with jobIndex =1 to numJobs will cover all scan iterations.
            %
            % See "BatchJobsExample.sh" for an example implementation.
            %
            % Inputs:
            % * qkdInput: A full QKDSolverInput with all required modules
            %   passed to MainIteration.
            % * basePath: string for the absolute or relative path to the
            %   folder where result files should be saved. If you construct
            %   this path in Matlab, we recommend using the function
            %   FULLFILE for OS agnostic file paths.
            % * baseFileName: The base file name used for output. The file
            %   will be named 'baseFileName_jobIndex.mat' to identify it.
            % * jobIndex: Positive 64 bit unsigned integer to distinguish
            %   this job from others. must be in the range 1 to
            %   numJobs.
            % * numJobs: positive integer for the total number of jobs
            %   the problem is split over. Must be smaller than or equal to
            %   the number of points to scan over.
            %
            % See also: MainIteration, QKDSolverInput, fullfile
            arguments
                qkdInput (1,1) QKDSolverInput
                basePath (1,1) string {mustBeFolder}
                baseFileName (1,1) string {checkBaseFileName}
                jobIndex (1,1) uint64 {mustBePositive}
                numJobs (1,1) uint64 {mustBeGreaterThanOrEqual(numJobs,jobIndex),...
                    numJobsMustNotExceedNumScanParams(numJobs,qkdInput)}
            end

            % Divide up the points to calculate for each job. If
            % numJobs doesn't evenly divide qkdInput.totalIterations,
            % then the first 'remainder' jobs each get an extra point to
            % calculate.
            numIterations = qkdInput.totalIterations;

            basePointsPerJob = idivide(numIterations,numJobs,"floor");
            remainder = numIterations - basePointsPerJob*numJobs;

            if jobIndex <= remainder
                % the first jobs get the extra remainder points
                startIndex = (basePointsPerJob + 1) * (jobIndex-1) + 1;
                endIndex   = startIndex + basePointsPerJob;
            else
                startIndex = basePointsPerJob * (jobIndex -1) + remainder + 1;
                endIndex = startIndex + basePointsPerJob - 1;
            end

            JobConstructor.runJob(qkdInput,basePath,baseFileName,jobIndex,startIndex:endIndex);
        end



        %% Compiling Jobs
        function compileJobs(basePath,baseFileName,clearJobFiles)
            % Compiles the files from individual job runs into a single
            % output file. COMPILEJOBS looks for all files that match the
            % naming scheme "baseFileName_#.mat" in the specified
            % directory. If no files are found, then a "NoMatchingFiles"
            % error is thrown.
            % 
            % Each, file must contain the variables:
            % * qkdInput: The QKDSolverInput used for all the runs. This is
            %   assumed to be the same across all files.
            % * results: The n x 1 struct array containing the results from
            %   each job. (It's the struct returned from calling
            %   MainIteration).
            % * scanIndexList: The n x 1 array that associates each entry
            %   in the results array with it's total scan index from
            %   qkdInput.
            %   
            % Compile jobs will concatenate all the results  and
            % scanIndexList arrays then sort them in ascending order using
            % scanIndexList.
            %
            % Inputs:
            % * basePath: String for the absolute or relative path to the
            %   folder where result files are saved. If you construct
            %   this path in Matlab, we recommend using the function
            %   FULLFILE for OS agnostic file paths.
            % * baseFileName: The base file name used for output. All files
            %   named 'baseFileName_*.mat' will be included.
            % * clearJobFiles (false): scalar logical. If true, the
            %   individual job files will be deleted after the final
            %   compiled file is created. If there are too few, too many,
            %   or duplicate entries in the concatenated scanIndexList,
            %   then a "CompilingIndexingProblem" warning is given and
            %   compileJobs will refuse to delete the individual job files.
            % See also: fullfile
            arguments
                basePath (1,1) string {mustBeFolder}
                baseFileName (1,1) string {checkBaseFileName}
                clearJobFiles (1,1) logical = false;
            end

            %% get and load all files that match the specified pattern.
            fileNamePattern = sprintf("%s_*.mat",baseFileName);
            files = dir(fullfile(basePath,fileNamePattern));

            % remove any directories that might have sneaked in.
            files = files(~[files(:).isdir]);
            % filter to 'fileName_#.mat'
            files = files(matches({files(:).name},...
                regexpPattern(regexptranslate("escape",baseFileName)+"_\d*\.mat")));
            

            if isempty(files)
                throw(MException("JobConstructor:NoMatchingFiles",...
                    "No files match the provided path and file name convention."))
            end

            filePathAndNames = string(fullfile({files(:).folder},{files(:).name}));


            fileContents = struct("results",cell(numel(files),1), ...
                "qkdInput",cell(numel(files),1), ...
                "scanIndexList",cell(numel(files),1));

            for index = 1:numel(files)
                [tmpContents.results,...
                    tmpContents.qkdInput,...
                    tmpContents.scanIndexList] ...
                    = checkJobFileFormat(load(filePathAndNames(index)));

                fileContents(index) = tmpContents;
            end

            %% compile results together
            % Ordered by the complete scanIndexList.
            qkdInput = fileContents(1).qkdInput;

            results = vertcat(fileContents(:).results);
            scanIndexList = vertcat(fileContents(:).scanIndexList);

            % sort the results by scanIndexList
            [sortedIndexes,sortOrder] = sort(scanIndexList);

            results = results(sortOrder);


            finalFile = fullfile(basePath,sprintf("%s.mat",baseFileName));
            save(finalFile,"results","qkdInput");


            % check to make sure the indexes are the indexes from 1 to the
            % total number of iterations. If not, we give a warning and
            % refuse to delete the individual job files.

            correctCompiledIndexing = isequal(sortedIndexes,(1:qkdInput.totalIterations)');
            if ~correctCompiledIndexing
                warning("JobConstructor:CompilingIndexingProblem",...
                    "Unexpected behaviour while compiling results. " + ...
                    "The compiled results from the job files do not " + ...
                    "match the indexing 1:qkdInput.totalIterations. " + ...
                    "The individual job files will not be deleted.");
            end


            %% remove all job files if the users asks to.
            if clearJobFiles && correctCompiledIndexing
                filePathAndNames = convertStringsToChars(filePathAndNames);
                if isa(filePathAndNames,"cell")
                    delete(filePathAndNames{:})
                else
                    delete(filePathAndNames)
                end
            end
        end
    end
end

%% validation functions

function mustBeInTotalScanParamRange(indexes,qkdInput)
if any(indexes>qkdInput.totalIterations)
    throwAsCaller(MException("JobConstructor:IndexOutOfScanRange",...
        "An index exceeds the total number of interations required to " + ...
        "enumerate all scan parameter combinations."))
end
end


function numJobsMustNotExceedNumScanParams(numJobs,qkdInput)
if numJobs > qkdInput.totalIterations
    throwAsCaller(MException("JobConstructor:NumJobsExceedsTotalScanIterations",...
        "The total number of Jobs exceeds the qkdInput's total number of " + ...
        "iterations to scan over."))
end
end

function [results,qkdInput,scanIndexList] = checkJobFileFormat(fileStruct)
arguments (Input)
    fileStruct (1,1) struct {mustHaveFields(fileStruct, ["results","qkdInput","scanIndexList"])}
end
arguments (Output)
    results (:,1) struct {mustHaveFields(results,["debugInfo","keyRate","currentParams"]),mustBeNonempty}
    qkdInput (1,1) QKDSolverInput
    scanIndexList (:,1) uint64 {mustBeNonzero,mustBeNonempty} % equal size at end of func.
end

results = fileStruct.results;
qkdInput = fileStruct.qkdInput;
scanIndexList = fileStruct.scanIndexList;

% Can't place in output arguments block. So here it goes
mustBeEqualSize(scanIndexList(:),results(:))

end

function checkBaseFileName(baseFileName)
% check for the base file name to ensure that it doesn't end with '_#'.
if any(endsWith(baseFileName,regexpPattern("_\d+")),"all")
    throwAsCaller(MException("JobConstructor","InvalidFileName",...
        "Base file names cannot end with '_#', for any number."))
end
end