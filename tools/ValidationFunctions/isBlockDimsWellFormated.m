function val = isBlockDimsWellFormated(blockDim)

%Case 1, Block dim is the default value (scalar and nan)
if isequaln(blockDim,nan)
    val = true;
    return
end

% Case 2, block dim is a vector of positive integers
if any(isnan(blockDim),"all")
    val = false;
    return
end
val = isvector(blockDim) && all(isfinite(blockDim))... %vector and finite
    && isequal(fix(blockDim),blockDim)... %integer
    && all(blockDim>0); % non negative

end