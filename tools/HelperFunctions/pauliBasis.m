function unitaryTrans = pauliBasis(basis,revY)
% pauliBasis returns the unitary 2 x 2 matrix that transforms the eigen
% states of the Pauli Z operator (H/1,V/2) to the eigen states of the Pauli
% X (D/+,A/-) or Y (R,L) operators.
%
% Note: Unlike, the functions present in Rotations, Qudit, and Coherent
% which are based on rotations around the Bloch sphere, pauliBasis does not
% pick up any additional phase factors.
%
% Inputs:
% * basis: Integer which selects which basis the unitary should transform
%   the Z eigen basis to: 1=Z (H/1,V/2) (identity), 2=X (D/+,A/-), 3=Y (R,L).
% * revY (false): Logical to reverse the assignment for the Pauli Y eigen
%   vectors. In other words, revY = true: U|1> = |R>, U|2> = |L>, revY =
%   false: U|1> = |L>, U|2> = |R>. Has no effect for other basis choices.
%
% See also: Rotations, Qudit, Coherent
arguments
    basis (1,1) double {mustBeMember(basis,[1,2,3])}
    revY (1,1) logical  = false
end
switch basis
    case 1
        %Z basis, no transformation needed
        unitaryTrans = eye(2);
    case 2
        %X basis, apply Hadamard
        unitaryTrans = [1,1;1,-1]/sqrt(2); %U|1> = |+>, U|2> = |->
    case 3
        % Y basis
        if revY
            unitaryTrans = [1,1;-1i,1i]/sqrt(2);% U|1> = |L>, U|2> = |R>
        else
            unitaryTrans = [1,1;1i,-1i]/sqrt(2);% U|1> = |R>, U|2> = |L>
        end
end
end