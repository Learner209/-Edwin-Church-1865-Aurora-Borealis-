MAGNETIC_MOMENT = 8 * 10 ^ 22; % magnetic moment of Earth
RE = 6371 * 1000; % Earth radius
MU = 4 * pi * 10 ^ (-7); % vaccum pemeatability
B0 = 3.07 * 10 ^ (-5); % Magnetic field above the equator
BASICE = 1.602176565 * 10 ^ (-19); % Basic e
CVELOCITY = 3 * 10 ^ 8 % c Velocity
MSTILL = 1.672621777 * 10 ^ (-27); % still Mass

centerX = 0;
centerY = 0;
centerZ = 0;
radius = RE;
[ballX, ballY, ballZ]=sphere(50);
ballX = centerX + radius * ballX;    
ballY = centerY + radius * ballY;
ballZ = centerZ + radius * ballZ;


cnt = 200;
streamlineCnt = 120;
scale = 50;
x =linspace(-scale * RE, 10 * scale * RE, cnt);
y =linspace(-scale * RE, scale * RE, cnt);
z =linspace(-scale * RE, scale * RE, cnt);

beta = pi / 2 * 0.9;
thetacnt = 6;
theta = linspace(- 1 * pi / 6 * 0.8, 1 * pi / 6 * 1.2, streamlineCnt);
Ax = RE * cos(beta) * cos(theta);
Ay = RE * cos(beta) * sin(theta);
Az = linspace(RE * sin(beta), RE * sin(beta), streamlineCnt);
startX = Ax;
startY = Ay;
startZ = Az;


[X, Y, Z] = meshgrid(x, y, z);
re = (X .^ 2 + Y .^ 2 + Z .^ 2) .^ (1 / 2);
% Sun exist beyond -20Re on the x_axis
rs = ((X + ones(size(X)) * 20 * RE).^ 2 + Y.^ 2 + Z.^ 2) .^ (1 / 2);

U = zeros(size(X));
V = zeros(size(Y));
W = zeros(size(Z));
% The dipole magnetic field of earth
U = U - MAGNETIC_MOMENT * MU / (4 * pi * re .^ 3) .* Z .* X * 3;
V = V - MAGNETIC_MOMENT * MU / (4 * pi * re .^ 3) .* Z .* Y * 3;
W = W - MAGNETIC_MOMENT * MU / (4 * pi * re .^ 3) .*(2 * Z .^ 2 - X .^ 2 - Y .^ 2);
% The dipole magnetic field of solar wind
U = U - 1 * MAGNETIC_MOMENT * MU / (4 * pi * rs .^ 3) .* Z .* (X + ones(size(X)) * 20 * RE) * 3;
V = V - 1 * MAGNETIC_MOMENT * MU / (4 * pi * rs .^ 3) .* Z .* Y * 3;
W = W - 1 * MAGNETIC_MOMENT * MU / (4 * pi * rs .^ 3) .*(2 * Z .^ 2 - (X + ones(size(X)) * 20 * RE) .^ 2 - Y .^ 2);

figure
surf(ballX, ballY, ballZ);

verts2 = stream3(X, Y, Z, U, V, W, startX, startY, -startZ);
lineobj = streamline(verts2);

view(3);

%Set initialized values
ENERGY = 10 * 10 ^ 6;

theta = pi / 6;
beta = pi / 4;
X0 = RE * cos(beta / 2) * sin(theta * 2);
Y0 = RE * cos(beta / 2) * cos(theta * 2);
Z0 = RE * sin(beta / 2);
X1 = -3 * 10000 * 0.995 * cos(beta) * sin(theta);
Y1 = 3 * 10000 * 0.995 * cos(beta) * cos(theta);
Z1 = 3 * 10000 * 0.995 * sin(beta);
[x, y] = ode45(@f, [0 1 * 10 ^ (-13)], [X0, Y0, Z0, X1, Y1, Z1]);

figure
surf(ballX, ballY, ballZ);
%disp(y(:,1));
plot3(y(:,1), y(:,2), y(:,3));
axis("tight");

function dy = f(x,y) % y contains m and q, 11 variables in total
    CVELOCITY = 3 * 10 ^ 8; %velocity of light
    Me = 1.672621777 * 10 ^ (-27); % still Mass
    BASICE = 1.602176565 * 10 ^ (-19); % Basic e
    
    %[azimuth,elevation,r] = cart2sph(y(1), y(2), y(3));
    %[XYZ,H,D,I,F] = wrldmagm(min(r, 85000), elevation, azimuth, decyear(2020,7,4),'Custom','WMM2020.COF');     
    B = magneticStrength(y(1), y(2), y(3));
%     B = zeros(3, 1);
%     B(1) = XYZ(1);
%     B(2) = XYZ(2);
%     B(3) = XYZ(3);

    dy = zeros(6, 1);
    dy(1) = y(4);
    dy(2) = y(5);
    dy(3) = y(6);
    dy(4) = BASICE / Me * (1 - (y(4) ^ 2 + y(5) ^ 2 + y(6) ^ 2) / CVELOCITY ^ 2) ^ (1/2) * (y(5) * B(3) - y(6) * B(2));
    dy(5) = BASICE / Me * (1 - (y(4) ^ 2 + y(5) ^ 2 + y(6) ^ 2) / CVELOCITY ^ 2) ^ (1/2) * (y(6) * B(1) - y(4) * B(3));
    dy(6) = BASICE / Me * (1 - (y(4) ^ 2 + y(5) ^ 2 + y(6) ^ 2) / CVELOCITY ^ 2) ^ (1/2) * (y(4) * B(2) - y(5) * B(1));
    disp([x]);
end

function B = magneticForce(x, y, z, magnetic_moment, mu, RE)

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
    U = U - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .* Z .* (X + ones(size(X)) * 20 * RE) * 3;
    V = V - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .* Z .* Y * 3;
    W = W - 1 * magnetic_moment * mu / (4 * pi * rs .^ 3) .*(2 * Z .^ 2 - (X + ones(size(X)) * 20 * RE) .^ 2 - Y .^ 2);
    B = zeros(1, 3);
    B(1, 1) = U;
    B(1, 2) = V;
    B(1, 3) = W;
end
