function value = isStructOf(inputStruct,type)
% isStructOf Validation function to determine if the given struct is of a
% specific type.
%
% Inputs:
% * inputStruct: The struct to test
% * type: The type to compare the struct against
%
% See also isCellOf
names = fieldnames(inputStruct);
value = ~any(cellfun(@(x) ~isa(inputStruct.(x),type),names));
end