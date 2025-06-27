%% Algorithm 5.4
clear
close all
clc

%% All dimensions in km, angles in deg, and angular rates in rad/s
mu = 398600;
f = 1/298.257223563; % Flattening of the Ellipsoid
OMEGA_E = 7.292115e-5; % Earth's Rotation Rate [rad/s]
OMEGA_vec = [0;0;OMEGA_E]; % [rad/s]
Re = 6378.1370; % [Km] Earth's Semi Major Axis, sometimes called a

%% Observer Location
theta = 300; % [deg] local sidereal time, not true anomaly!
phi = 60; % [deg] Latitude
H = 0; % [Km] Altitude

%% Satellite Measurements
rho = 2551; % [Km] Range
A = 90; % [deg] Azimuth
a = 30; % [deg] Elevation
rho_d = 0; % [Km/s] Range Rate
A_d = 1.973e-3; % [rad/s] Azimuth Rate
a_d = 9.864e-4; % [rad/s] Elevation Rate

%% Algorithm 5.4
R_vec = (Re/sqrt(1-(2*f-f^2)*sind(phi)^2)+H) * [cosd(phi)*cosd(theta); cosd(phi)*sind(theta); (1-f)^2*sind(phi)];
delta = asind(cosd(phi)*cosd(A)*cosd(a) + sind(phi)*sind(a));
if(A>0 & A<180)
    h = 360 - acosd((cosd(phi)*sind(a)-sind(phi)*cosd(A)*cosd(a))/cosd(delta));
else
    h = acosd((cosd(phi)*sind(a)-sind(phi)*cosd(A)*cosd(a))/cosd(delta));
end
alpha = theta - h;
rho_hat_vec = [cosd(delta)*cosd(alpha); cosd(delta)*sind(alpha); sind(delta)];
r_vec = R_vec + rho*rho_hat_vec;
R_d_vec = cross(OMEGA_vec, R_vec);
delta_d = 1/cosd(delta)*(-A_d*cosd(phi)*sind(A)*cosd(a) + a_d*(sind(phi)*cosd(a)-cosd(phi)*cosd(A)*sind(a)));
alpha_d = OMEGA_E + (A_d*cosd(A)*cosd(a)-a_d*sind(A)*sind(a)+delta_d*sind(A)*cosd(a)*tand(delta)) / (cosd(phi)*sind(a)-sind(phi)*cosd(A)*cosd(a));
rho_hat_d_vec = [-alpha_d*sind(alpha)*cosd(delta)-delta_d*cosd(alpha)*sind(delta);
                 alpha_d*cosd(alpha)*cosd(delta)-delta_d*sind(alpha)*sind(delta);
                 delta_d*cosd(delta)];
v_vec = R_d_vec + rho_d*rho_hat_vec + rho*rho_hat_d_vec;

%% You can now use r & v to proceed with algorithm 4.2
[a,e,i,OMEGA,omega,theta] = orbital_elements(r_vec,v_vec)
