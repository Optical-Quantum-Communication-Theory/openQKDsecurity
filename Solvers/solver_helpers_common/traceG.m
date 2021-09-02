function fval = traceG(rho,krausOperators)

    % for the case there is a post-selection map.

    gRho = krausFunc(rho,krausOperators); % calculate G(\rho).
    dim = size(gRho,1); % get the dimension of G(\rho).
    eigMin = lambda_min(gRho); % check the minimum eigenvalue of this density matrix
    if eigMin <= 0
       % if the eigenvalue is less than one, do a perturbation by using a
        % pinching channel.
       epsilon = (1e-14-eigMin)*dim;
       gRho = (1-epsilon)*gRho + epsilon*eye(dim)/dim;
    end

    fval = real(trace(gRho)); % calculate the gain

end