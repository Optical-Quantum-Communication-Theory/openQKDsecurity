%% FUNCTION NAME: removeLinearDependence
% Removes linearly dependent matrices from a set of square matrices.

function [independentSet,independentCols] = removeLinearDependence(dependentSet)

    dim = size(dependentSet{1},1);

    % Reshape the set into a matrix
    reshapedDependentSet = cell(size(dependentSet));
    for iEntry = 1:numel(dependentSet)
        reshapedDependentSet{iEntry} = reshape(dependentSet{iEntry},1,dim*dim);
    end
    reshapedDependentSet = transpose(cell2mat(reshapedDependentSet));

    % Remove linearly dependent matrices
    reshaped_set=size(reshapedDependentSet);
    [~,independentCols] = rref(reshapedDependentSet,10^-20);
%     isequal(independentCols,independentCols1)
    independentSet = dependentSet(independentCols);

end