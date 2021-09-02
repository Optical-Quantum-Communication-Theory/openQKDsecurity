function zRho = zMap(rho,keyMap)
    zRho = 0;
    for i = 1 : numel(keyMap)
        zRho = zRho + keyMap{i}*rho*keyMap{i};
    end
end