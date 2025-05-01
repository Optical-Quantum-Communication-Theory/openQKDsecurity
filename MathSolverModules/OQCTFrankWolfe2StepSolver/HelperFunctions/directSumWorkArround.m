function rhoCVX = directSumWorkArround(rhoCVX,rhoBlocks)
% CVX's has a bug with blkdiag, so we had to write our own work around.
% This manually builds the direct sum by placing the blocks on the
% diagonal.
%
% Works only for square n x n matrices!
%
% Inputs:
% * rhoCVX: n x n matrix defined by CVX
% * rhoBlocks: cell array containing the individual block matrices, such
%   that directsum_i rhoBlocks{i} = rhoCVX 
%
% Outputs:
% * rhoCVX: n x n matrix defined by CVX after performing direct sum

%Define stride length
runningSum = 0;

%Perform direct sum by running over block-sizes
for index = 1:numel(rhoBlocks)
    %Find size of current block
    thisBlockSize = size(rhoBlocks{index},1);
    %Assigns entries of "big" rho defined by CVX to equal corresponding
    %blocks
    rhoCVX(runningSum+1:runningSum+thisBlockSize,runningSum+1:runningSum+thisBlockSize) = rhoBlocks{index};
    %Update stride length
    runningSum = runningSum + size(rhoBlocks{index},1);
end
end