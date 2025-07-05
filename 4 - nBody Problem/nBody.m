clear
clc
close all

%% Problem Data
G = 6.6743e-11; % [m^3/kg/s^2]
mu = G*1.989e30; % assume mu depends only on Sun

tfinal = 365*24*3600; % 1*12*30*24*3600; % 67000;
dt = 30;
t = 0:dt:tfinal;

%% Bodies Data
obj1.filename = "./Images/sun.png";
obj1.radius = 696e6 * 30;
r1_0 = [0; 0; 0];
v1_0 = [0e3; 0e3; 0e3];
m1 = 1.989e30; % [Kg]

obj2.filename = "./Images/mercury.png";
obj2.radius = 2440e3 * 900;
a_mercury = 57.91e9; % [m] semi-major axis
e_mercury = 0.2056; % eccentricity
i_mercury = 7.00 *pi/180; % Inclination of Orbit to the Ecliptic Plane
r_mercury = a_mercury*(1+e_mercury);
v_mercury = sqrt(mu*(2/r_mercury-1/a_mercury));
r2_0 = [-r_mercury; 0; 0]; % apogee
v2_0 = [0; -v_mercury; 0];
r2_0 = eul2rotm([0,i_mercury,0]) * r2_0;
v2_0 = eul2rotm([0,i_mercury,0]) * v2_0;
m2 = 330.2e21; % [Kg]

obj3.filename = "./Images/venus.jpg";
obj3.radius = 6052e3 * 900;
a_venus = 108.2e9; % [m] semi-major axis
e_venus = 0.0067; % eccentricity
i_venus = 3.39 *pi/180; % Inclination of Orbit to the Ecliptic Plane
r_venus = a_venus*(1+e_venus); % apogee
v_venus = sqrt(mu*(2/r_venus-1/a_venus));
r3_0 = [-r_venus; 0; 0];
v3_0 = [0; -v_venus; 0];
r3_0 = eul2rotm([0,i_venus,0]) * r3_0;
v3_0 = eul2rotm([0,i_venus,0]) * v3_0;
m3 = 4.869e24; % [Kg]

obj4.filename = "./Images/earth.png";
obj4.radius = 6378e3 * 900;
a_earth = 149.6e9; % [m] semi-major axis
e_earth = 0.0167; % eccentricity
i_earth = 0.00 *pi/180; % Inclination of Orbit to the Ecliptic Plane
r_earth = a_earth*(1+e_earth); % apogee
v_earth = sqrt(mu*(2/r_earth-1/a_earth));
r4_0 = [-r_earth; 0; 0];
v4_0 = [0; -v_earth; 0];
r4_0 = eul2rotm([0,i_earth,0]) * r4_0;
v4_0 = eul2rotm([0,i_earth,0]) * v4_0;
m4 = 5.974e24; % [Kg]

obj5.filename = "./Images/mars.jpg";
obj5.radius = 3396e3 * 900;
a_mars = 227.9e9; % [m] semi-major axis
e_mars = 0.0935; % eccentricity
r_mars = a_mars*(1+e_mars); % apogee
v_mars = sqrt(mu*(2/r_mars-1/a_mars));
r5_0 = [-r_mars; 0; 0];
v5_0 = [0; -v_mars; 0];
m5 = 641.9e21; % [Kg]

obj_vec = [obj1, obj2, obj3, obj4, obj5];
states_0 = [r1_0; r2_0; r3_0; r4_0; r5_0; v1_0; v2_0; v3_0; v4_0; v5_0];
mass_vec = [m1, m2, m3, m4, m5];

%% Solve
tic;
[t, states] = ode45(@(t, states_0) fdot(t, states_0, mass_vec), t, states_0);
toc;

N = length(mass_vec);

R = states(:, 1:3*N);
% RG = [sum(mass_vec.*R(:,1:3:3*N), 2), sum(mass_vec.*R(:,2:3:3*N), 2), sum(mass_vec.*R(:,3:3:3*N), 2)]  /  sum(mass_vec); % Position of Centre of Gravity
% R = [R, RG];
clear RG

V = [];
% V = states(:,3*N+1:2*3*N);
% VG = [sum(mass_vec.*V(:,1:3:3*N), 2), sum(mass_vec.*V(:,2:3:3*N), 2), sum(mass_vec.*V(:,3:3:3*N), 2)]  /  sum(mass_vec); % Velocity of Centre of Gravity
% V = [V, VG];
% clear VG

clear states

A = [];
% A = nan(size(V));
% for i = 1:N
%     array = 1:N; % indices of all bodies except current one
%     array = array(array~=i);
%     R_i = R(:,i);
%     A(:,N+i) = sum( G*mass_vec(array).*(R(:,array) - R_i) ./ vecnorm(R(:,array) - R_i).^3, 2);
% end

R = R.';
V = V.';
A = A.';

%% Animate & Plot
fig1 = figure('Position', [1, 1, 1366, 728]);
ax1 = gca;
axis equal
% pbaspect([1,1,1]);
xlim([-2.5e11, 2.5e11]);
ylim([-2.4e11, 2.4e11]);
zlim([-7e10, 7e10]);
view([-50, 15]);
hold on; grid on;
xlabel("x [Km]", 'Interpreter', 'latex'); ylabel("y [Km]", 'Interpreter', 'latex'); zlabel("z [Km]", 'Interpreter', 'latex');
title("The n-Body Problem", 'Interpreter', 'latex');
simulate(t, 86400/dt, R, V, A, ["r", "g", "m", "b", "c"], ["Sun", "Mercury", "venus", "Earth", "Mars"], ax1, obj_vec);
% max(vecnorm(V(1:3,:),2,1))
