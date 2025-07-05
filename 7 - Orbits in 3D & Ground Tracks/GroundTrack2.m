clear
clc
close all

%% Earth Data
R_earth = 6378; % [km]
mu = 398600; % [km^3/s^2]
J2 = 1.08263e-03;

%% Problem Data
rp = 6700; % [km]
ra = 10000; % [km]
theta0 = 230; % [deg]
OMEGA0 = 270; % [deg]
i0 = 60; % [deg]
omega0 = 45; % [deg]
delta_t = 45*60; % [sec]

a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
T = 2*pi/sqrt(mu)*a^(3/2)
OMEGA_d = -3/2*sqrt(mu)*J2*R_earth^2/((1-e^2)^2*a^3.5)*cosd(i0) *180/pi % [deg/s]
omega_d = -3/2*sqrt(mu)*J2*R_earth^2/((1-e^2)^2*a^3.5)*(sind(i0)^2-2) *180/pi % [deg/s]

E0 = 2*atan(sqrt((1-e)/(1+e))*tand(theta0/2)) % [rad]
M0 = E0 - e*sin(E0) % [rad]
t0 = M0*T/2/pi

%% at time t0+delta_t
t = t0+delta_t
M = t*2*pi/T
E = newton_raphson(@(E) E-e*sin(E)-M, @(E) 1-e*cos(E), E0)
theta = 2*atand(sqrt((1+e)/(1-e))*tan(E/2)) % [deg]
OMEGA = OMEGA0 + OMEGA_d*delta_t
omega = omega0 + omega_d*delta_t
