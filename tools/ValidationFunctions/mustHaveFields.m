function mustHaveFields(inputStruct,fields)
% mustHaveField Throws an error if the struct (or struct array) does not
% contain all listed fields.
%
% Inputs:
% * inputStruct: The struct (or struct array) to check.
% * fields: String array of field names the struct must have.
fieldCheck = isfield(inputStruct,fields(:));
if ~all(fieldCheck,"all")
    throwAsCaller(MException("ValidationFunction:StructMissingFields",...
        "The struct is missing the field(s) '%s'.",strjoin(fields(~fieldCheck),"', '")))
end
end