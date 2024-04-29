function primalConstraints(cvxProb,rho,eqCon,ineqCon,vec1NormCon,mat1NormCon,linConTol)
%All of the primal problem constraints packaged into a single function
%call. Adds all of them to the given CVX problem using the newcnstr
%function. This function was split off from the FW2stepSolver for better
%readability.
primalLinearEqualityConstraints(cvxProb,rho,eqCon,linConTol);
primalLinearInequalityConstraints(cvxProb,rho,ineqCon,linConTol);
primalVectorOneNormConstraints(cvxProb,rho,vec1NormCon,linConTol);
primalMatrixOneNormConstraints(cvxProb,rho,mat1NormCon,linConTol);
end


function primalLinearEqualityConstraints(cvxProb,rho,equalityConstraints,linearConstraintTolerance)
%Adds the list of equality constraints using the newcnstr function from
%CVX.
for index = 1:numel(equalityConstraints)
    newcnstr(cvxProb,...
        real(trace(equalityConstraints(index).operator'*(rho))),...
        equalityConstraints(index).scalar + linearConstraintTolerance,...
        '<=');
    newcnstr(cvxProb,...
        real(trace(equalityConstraints(index).operator'*(rho))),...
        equalityConstraints(index).scalar - linearConstraintTolerance,...
        '>=');
end
end


function primalLinearInequalityConstraints(cvxProb,rho,inequalityConstraints,linearConstraintTolerance)
%Adds the list of inequality constraints using the newcnstr function from
%CVX.
for index = 1:numel(inequalityConstraints)
    if ~isinf(inequalityConstraints(index).upperBound)
        newcnstr(cvxProb,...
            real(trace(inequalityConstraints(index).operator'*(rho))),...
            inequalityConstraints(index).upperBound + linearConstraintTolerance,...
            '<=');
    end
    if ~isinf(inequalityConstraints(index).lowerBound)
        newcnstr(cvxProb,...
            real(trace(inequalityConstraints(index).operator'*(rho))),...
            inequalityConstraints(index).lowerBound - linearConstraintTolerance,...
            '>=');
    end
end
end


function primalVectorOneNormConstraints(cvxProb,rho,vec1NormCons,linConTol)
%Adds the list of vector one norm constraints using the one norm and the
%newcnstr function from CVX.
for index = 1 : numel(vec1NormCons)
    newcnstr(cvxProb,...
        norm(vecstatmap(rho, vec1NormCons(index).operators)...
        - vec1NormCons(index).vector,1),... 
        vec1NormCons(index).scalar + linConTol,...
        '<=');
end
end


function primalMatrixOneNormConstraints(cvxProb,rho,mat1NormCons,linConTol)
%Adds the list of matrix one norm constraints using the nuclear norm (trace
%norm, i.e. those of the form ||map(rho) - matrix||_1 <= scalar) and the
%newcnstr function from CVX.
for index = 1 : numel(mat1NormCons)
    mappedrho = ApplyMap(rho, mat1NormCons(index).choiMatrix);
    newcnstr(cvxProb,...
        norm_nuc(mappedrho - mat1NormCons(index).operator),...
        mat1NormCons(index).scalar+linConTol,...
        '<=');
end
end


%% helper functions

function distribution = vecstatmap(rho, observables)
% Function that maps expectations onto the diagonal entries of a matrix.
% Used for one-norm constraints.
for index = numel(observables):-1:1
    distribution(index) = trace(observables{index}'*rho);
end
% matlab constructs it in the wrong shape when using reverse indexing. We
% just fix that rate here. (Transpose, NOT conjugate).
distribution = distribution.'; 
end