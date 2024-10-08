clear
% A little script to check and see if all the components are installed and
% added to the path by calling function from each.

% Matlab is obvious so we don't need to try that one.

numProblems =0;

%% CVX
disp("Checking if CVX is installed")
try
    %very simple CVX problem
    cvx_begin quiet
    variables x(2)
    minimize norm(x)
    subject to
    x(1) +x(2) <= 5;
    x(1) -x(2) >= 2;
    cvx_end
    disp("Passed")
catch err
    numProblems = numProblems+1;
    if ~checkCantFindFunction(err.identifier,"CVX")
        disp("Some VERY unexpected error occured. Here's the stack trace.")
        rethrow(err);
    end
end

%% Qetlab
disp("Checking if Qetlab is installed")
QetlabFound = false;
try
    iden(1,false);
    QetlabFound = true;
    disp("Passed")
catch err
    numProblems = numProblems+1;
    if ~checkCantFindFunction(err.identifier,"Qetlab")
        disp("Some VERY unexpected error occured. Here's the stack trace.")
        rethrow(err);
    end
end

if QetlabFound
    % Qetlab version 0.9 and lower has major bugs with using choi matrices.
    % Unfortunetly, Qetlab has not updated their stable branch on their website
    % so we have to direct people to their github.
    disp("Checking if Qetlab version is above 0.9")
    % The choi matrix we use for the test comes from this channel on 2x2 rho.
    % testChannel = @(rho) blkdiag(0.25*rho,0.75*trace(rho));
    testChoiMat = zeros(6,6);
    testChoiMat(1,1) = 0.25;
    testChoiMat(1,5) = 0.25;
    testChoiMat(5,1) = 0.25;
    testChoiMat(5,5) = 0.25;
    testChoiMat(3,3) = 0.75;
    testChoiMat(6,6) = 0.75;

    qetLabChoiMat =  PartialMap(MaxEntangled(2,false,false)*MaxEntangled(2,false,false)',...
        testChoiMat,2,[2,2]);

    % For Qetlab> 0.9, these should be EXACTLY the same.
    if ~all(testChoiMat == qetLabChoiMat,"all")
        numProblems = numProblems+1;
        warning("Qetlab's PartialMap function failed for our test Choi matrix." + ...
            " This is likely because you are using Qetlab version 0.9 or lower." + ...
            " Please go to https://github.com/nathanieljohnston/QETLAB for the" + ...
            " latest version of Qetlab.")
    else
        disp("Passed")
    end

end

%% OpenQKDSecurity
disp("Checking if OpenQKDSecurity is installed")
try
    debugInfo = DebugInfo();
    disp("Passed")
catch err
    numProblems = numProblems+1;
    if ~checkCantFindFunction(err.identifier,"OpenQKDSecurity")
        disp("Some VERY unexpected error occured. Here's the stack trace.")
        rethrow(err);
    end
end

%% All done!

if numProblems == 0
    disp("Completed with no errors.")
else
    fprintf("Completed with %d problem(s) to address.\n",numProblems)
end

clear

function isTheError = checkCantFindFunction(identifier,packageName)
isTheError = false;
if isequal(identifier, "MATLAB:UndefinedFunction")
    isTheError = true;
    warning("It seems that the package %s is missing. Plsease ensure" + ...
        " that the package is installed and all folders and subfolders" + ...
        " are on the path. If you encountered this error after restarting" + ...
        " Matlab, then you likely didn't save your path setup and you" + ...
        " have to redo them.\n",packageName);
end
end