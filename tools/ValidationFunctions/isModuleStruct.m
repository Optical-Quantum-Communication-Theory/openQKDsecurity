function value = isModuleStruct(module)
% isModuleStruct Validation function to determine if the given struct is a
% module struct by testing if it has the field "moduleFunc".
%
% Inputs:
% module: The struct to test
if ~isstruct(module)
    value = false;
    return
end
value = isfield(module,"moduleFunc");
end