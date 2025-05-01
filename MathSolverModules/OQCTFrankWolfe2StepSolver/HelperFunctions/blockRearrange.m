function [P, newDims] = blockRearrange(a, b)
% blockRearrange Gives the matrix required to re-arrange a tensor product 
% of block diagonal matrices into its best block diagonal form.
% returns a matrix P such that, if M is a tensor product of two block
% diagonal matrices A and B, whose block sizes are represented by the 
% vectors a and b, then P*M*P' will be a block diagonal matrix whose
% block sizes are represented by the vector kron(a,b).
% Written by Scott Johnstun
%
%
% Inputs:
% * a: The block dimensions of Alice
% * b: The block dimensions of Bob
%
% Outputs:
% * P: A matrix that gives the rules for how to rearrange rows to maximize
%   block diagonal structure. If M is a tensor product of block diagonal
%   matrices A and B, then P*M*P' will be equivalent to M, but with rows
%   and columns rearranged to preserve maximal block diagonal structure.
% * newDims: The block diagonal dimensions of the new matrix P*M*P'.
%   Computed as kron(a,b).
%
% Underlying function
% $n$ and $m$ are the dimensions of $\vec{a}$ and $\vec{b}$ respectively.
%
% $\text{SWAP}_{l,r} = \sum_{i=1}^l \sum_{j=1}^r (\vec{e}_j \otimes
% \vec{e}_i})(\vec{e}_i \otimes \vec{e}_j)^\tran$
%
% $P = \bigoplus_{i=1}^n \left(\bigoplus_{j=1}^m \text{SWAP}_{b_j,a_i}
% \right) \text{SWAP}_{a_i,\vec{1s}\cdot\vec{b}}$
arguments
    a (1,:) {mustBeInteger,mustBePositive,mustBeNonempty}
    b (1,:) {mustBeInteger,mustBePositive,mustBeNonempty}
end
% rewritten by John Burniston
dimA = sum(a,"all");
dimB = sum(b,"all");
dimP = dimA*dimB;

% Instead of concatenating a large number of blocks together, we start with
% small vectors that represent permutations, then at the end we expand that
% into a permutation matrix. It turns out, the easiest way to do this is to
% compute the inverse permutation vector,$\vec{pi}^{-1}$ then map it to the
% correct permutation matrix $P_\pi$ at the very end.
Y = cell(1,numel(a));
for indexA = 1:numel(a) %indexing over Alice's blocks
    if a(indexA) == 1
        % when a_i ==1, we don't have to permute Bob's systems for this
        % chunk's max block diagonal structure.
        Y{indexA} = 1:dimB;
    else
        % For Alice's indexAth block we compute the equivalent of
        % $\text{SWAP}^\dagger_{b_j,a_i}$ for $j=1,\dots,m$ and store them
        % in pi2.
        pi2inv = arrayfun(@(b_j)swapPermInv(b_j,a(indexA)),b,"UniformOutput",false);

        % If we we're computing  $bigoplus_{j=1}^m
        % \text{SWAP}^\dagger_{b_j,a_i}$  directly, then we could just call
        % blkdiag to group them together. However, using the permutation
        % vectors, we have to be a little more careful with our indices.
        % The permCat function will take care of this for us.
        pi2inv = permCat(pi2inv{:});

        % This outside loop is doing a similar thing for
        % $\text{SWAP}^\dagger_{a_i,\vec{1s}\cdot\vec{b}}$.

        % Let $\pi^{-1}_1$ =  swapPermInv(a(indexA),dimB). Because we're
        % computing $\pi^{-1}$, instead of storing $\pi_2 \circ \pi_1$, we
        % are storing $\pi^{-1}_1 \circ \pi^{-1}_2$.
        Y{indexA} = pi2inv(swapPermInv(a(indexA),dimB));
    end
    
end
% Just like before we concatenate the permutations using our special little
% function.
Y = permCat(Y{:});

% These next two lines map $\pi^{-1}$ to $P_\pi$.
P = speye(dimP);
P = P(Y,:)';

newDims = kron(a,b);

end


function perm = swapPermInv(m,n)
%computes the inverse permutation vector for swapping two systems.
perm = reshape(1:m*n,[m,n])';
perm = perm(:);
end

function perm = permCat(a)
arguments (Repeating)
    a (1,:) {mustBeNonempty}
end
shift=0;
numelA = cellfun(@numel,a);
perm = zeros(1,sum(numelA));
for index = 1:numel(a)
    perm(shift+1:shift+numelA(index)) = a{index}+shift;
    shift = shift + numelA(index);
end
end

%% Scott's original (though slightly buggy when the result should be identity)
% function [P, newDims] = blockRearrange(a, b)
% arguments
%     a (1,:) {mustBeInteger}
%     b (1,:) {mustBeInteger}
% end
% % reshape a and b to make sure they are row vectors
%     a = reshape(a, 1,[]);
%     b = reshape(b, 1,[]);
%     dimA = sum(a);
%     dimB = sum(b);
%     dimP = dimA * dimB;
%     P = zeros(dimP);
%     for i = 1 : length(a)
%         if a(i) == 1
%             P(sum(a(1:i-1))*dimB+1 : sum(a(1:i))*dimB, sum(a(1:i-1))*dimB+1 : sum(a(1:i))*dimB) = eye(dimB);
%         elseif a(i) > 1
%             Y = commutation(a(i), dimB);
%             P(sum(a(1:i-1))*dimB+1 : sum(a(1:i))*dimB, sum(a(1:i-1))*dimB+1 : sum(a(1:i))*dimB) = Y;
%         elseif a(i) < 1
%             warning("found an invalid dimension at position %d in Alice's array: %d", i, a(i))
%         end
%     end
%     newDims = kron(a,b);
% end
% 
% % returns the commutation matrix of size (m,n)
% % taken from
% % https://stackoverflow.com/questions/52569719/how-to-compute-commutation-matrix-in-matlab
% function Y = commutation(m,n)
%     I = reshape(1:m*n, [m, n]); % initialize a matrix of indices of size(A)
%     I = I'; % Transpose it
%     I = I(:); % vectorize the required indices
%     Y = speye(m*n); % Initialize an identity matrix
%     Y = Y(I,:); % Re-arrange the rows of the identity matrix
%     Y = Y';
% end