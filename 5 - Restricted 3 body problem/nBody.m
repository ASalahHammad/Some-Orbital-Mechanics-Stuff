clear all
clc
close all

%% Problem Data
G = 6.6743e-11; % [m^3/kg/s^2]
mu = G*(5.974e24 + 73.48e21);

tfinal = 1000; % 10*24*3600; % 67000;
dt = 1;
t = 0:dt:tfinal;

%% Bodies Data
obj1.filename = "./Images/earth.png";
obj1.radius = 6378e3 *20;
a_earth = 149.6e9; % [m] semi-major axis
e_earth = 0.0167; % eccentricity
i_earth = 0.00 *pi/180; % Inclination of Orbit to the Ecliptic Plane
r_earth = a_earth*(1+e_earth); % apogee
v_earth = sqrt(mu*(2/r_earth-1/a_earth));
r1_0 = [0; 0; 0];
v1_0 = [0; 0; 0];
r1_0 = eul2rotm([0,i_earth,0]) * r1_0;
v1_0 = eul2rotm([0,i_earth,0]) * v1_0;
m1 = 5.974e24; % [Kg]

obj2.filename = "./Images/moon.png";
obj2.radius = 1737e3 *10;
a_moon = 384.4e6; % [m] semi-major axis
e_moon = 0.0549; % eccentricity
i_moon = 5.145 *pi/180; % Inclination of Orbit to the Ecliptic Plane
r_moon = a_moon*(1+e_moon); % apogee
v_moon = sqrt(mu*(2/r_moon-1/a_moon));
r2_0 = [r_moon; 0; 0];
v2_0 = [0; v_moon; 0];
% r2_0 = eul2rotm([0,i_moon,0]) * r2_0;
% v2_0 = eul2rotm([0,i_moon,0]) * v2_0;
m2 = 73.48e21; % [Kg]

% obj3.filename = "./Images/moon.png";
% obj3.radius = 1737e3 *10;
% a_moon = 384.4e6; % [m] semi-major axis
% e_moon = 0.0549; % eccentricity
% i_moon = 5.145 *pi/180; % Inclination of Orbit to the Ecliptic Plane
% r_moon = a_moon*(1+e_moon); % apogee
% v_moon = sqrt(mu*(2/r_moon-1/a_moon));
% r3_0 = [3.264e8; 0e8; 0];
% v3_0 = [0; sqrt(G*m1/r3_0(1)); 0];
% % r2_0 = eul2rotm([0,i_moon,0]) * r2_0;
% % v2_0 = eul2rotm([0,i_moon,0]) * v2_0;
% m3 = 50; % [Kg]

%% Move Origin
% r1_0 = r1_0 + [2e8;0;0];
% r2_0 = r2_0 + [2e8;0;0];

obj_vec = [obj1, obj2];
states_0 = [r1_0; r2_0; v1_0; v2_0];
mass_vec = [m1, m2];

%% Solve
tic;
[t, states] = ode45(@(t, states_0) fdot(t, states_0, mass_vec), t, states_0);
toc;

N = length(mass_vec);

R = states(:, 1:3*N);
RG = [sum(mass_vec.*R(:,1:3:3*N), 2), sum(mass_vec.*R(:,2:3:3*N), 2), sum(mass_vec.*R(:,3:3:3*N), 2)]  /  sum(mass_vec); % Position of Centre of Gravity
R = [R, RG];
% clear RG

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
RG = RG.';
V = V.';
A = A.';

% Ro = (R(1:3,:) + R(4:6,:)) / 2;

%% Animate & Plot
fig1 = figure('Position', [1, 1, 1368, 728]);
ax1 = gca;
axis equal
% pbaspect([1,1,1]);
xlim([-4.6e8, 5.2e8]);
ylim([-4e8, 4e8]);
% zlim([-2e7, 4e7]);
view(2);
hold on; grid on;
xlabel("x [m]", 'Interpreter', 'latex'); ylabel("y [m]", 'Interpreter', 'latex'); zlabel("z [m]", 'Interpreter', 'latex');
title("The n-Body Problem", 'Interpreter', 'latex');
simulate(t, 5000, mass_vec, R(1:end-3,:), V, A, ["b", "r--"], ["Earth", "Moon"], ax1, obj_vec);
% max(vecnorm(V(1:3,:),2,1))
