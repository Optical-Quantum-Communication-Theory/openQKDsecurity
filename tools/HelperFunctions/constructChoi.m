function J = constructChoi(map, dimX,useSparseArrays)
% CONSTRUCTCHOI Constructs the Choi matrix representation of a  map from
% the spaces X to Y. Computes J(Phi) = sum_i,j |i><j| \otimes Phi(|i><j|).
%
% Inputs:
% * map: A function handle that applies Phi
% * dimX: The dimensions of the input space for the channel
% * useSparseArrays (false): Construct the Choi matrix using sparse arrays.
% 
% TODO:
% * Make an option to generate with sparse arrays
arguments
    map (1,1) function_handle
    dimX (1,1) double {mustBeInteger,mustBePositive}
    useSparseArrays (1,1) logical = false;
end

if useSparseArrays
    dimY = size(map(sparse(dimX,dimX)));
else
    dimY = size(map(zeros(dimX)));
end

% make sure the output is 2D and square or else this won't work.
if numel(dimY)>2 || dimY(1) ~= dimY(2)
    throw(MException("constructChoi:MapDimensionsNotSupported",...
        "The map must take square 2D matrices to square 2D matrices."))
end
dimY = dimY(1);

if useSparseArrays
    J = sparse(dimX*dimY,dimX*dimY);
    basisElement = sparse(dimX,dimX);
else
    J = zeros(dimX*dimY, dimX*dimY);
    basisElement = zeros(dimX,dimX);
end

% loop and create the Choi matrix
for iKet = 1 : dimX
    ketStart = dimY*(iKet-1)+1;
    KetStop = dimY*iKet;
    for jBra = 1 : dimX
        braStart = dimY*(jBra-1)+1;
        braStop = dimY*jBra;
        basisElement(iKet,jBra) = 1; % slightly faster than making it each time.
        J(ketStart:KetStop,braStart:braStop) =map(basisElement);
        basisElement(iKet,jBra) = 0;

        % functionally equivalent but takes way too much time to compute.
        % J = J + kron(basisElement, map(basisElement));
    end
end
end