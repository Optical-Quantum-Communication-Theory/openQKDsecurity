function value = isStructOf(inputStruct,type)
% isStructOf Validation function to determine if the given struct is of a
% specific type.
%
% Inputs:
% * inputStruct: The struct to test
% * type: The type to compare the struct against
%
% See also isCellOf
arguments
    inputStruct (1,1) struct
    type (:,1) string{mustBeNonempty,mustBeNonzeroLengthText}
end
value = ~any(structfun(@(x) ~isa(x,type),inputStruct),"all");
end