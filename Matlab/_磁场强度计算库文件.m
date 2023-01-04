function B = magneticStrength(x, y, z, magnetic_moment, mu, RE)

    if(~exist('mu','var'))
        mu = 4 * pi * 10 ^ (-7); % vaccum pemeatability;
    end
    if(~exist('Re','var'))
        RE = 6371 * 1000; % Earth radius
    end
    if(~exist('magnetic_moment','var'))
        magnetic_moment = 8 * 10 ^ 22; % magnetic moment of Earth
    end
    [X, Y, Z] = meshgrid(x, y, z);
    Re = (X .^ 2 + Y .^ 2 + Z .^ 2) .^ (1 / 2);
    % Sun exist beyond -20Re on the x_axis
    rs = ((X + ones(size(X)) * 20 * RE).^ 2 + Y.^ 2 + Z.^ 2) .^ (1 / 2);
    
    U = zeros(size(X));
    V = zeros(size(Y));
    W = zeros(size(Z));
    % The dipole magnetic field of earth
    U = U - magnetic_moment * mu / (4 * pi * Re .^ 3) .* Z .* X * 3;
    V = V - magnetic_moment * mu / (4 * pi * Re .^ 3) .* Z .* Y * 3;
    W = W - magnetic_moment * mu / (4 * pi * Re .^ 3) .*(2 * Z .^ 2 - X .^ 2 - Y .^ 2);
    % The dipole magnetic field of solar wind
    %U = U - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .* Z .* (X + ones(size(X)) * 20 * RE) * 3;
    %V = V - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .* Z .* Y * 3;
    %W = W - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .*(2 * Z .^ 2 - (X + ones(size(X)) * 20 * RE) .^ 2 - Y .^ 2);
    B = zeros(1, 3);
    B(1, 1) = U;
    B(1, 2) = V;
    B(1, 3) = W;
end
