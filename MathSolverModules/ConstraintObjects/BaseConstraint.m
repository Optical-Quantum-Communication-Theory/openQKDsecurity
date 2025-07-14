classdef (Abstract) BaseConstraint
    %BASECONSTRAINT Just a basic abstract class we can put universal
    %function handles for all our classes. Nothing really to see here.

    properties (Abstract = true, Dependent = true)
        rhoDim (1,1) double {mustBeNonnegative} % Size of the input system's hilbert space for this constraint.
    end
end

