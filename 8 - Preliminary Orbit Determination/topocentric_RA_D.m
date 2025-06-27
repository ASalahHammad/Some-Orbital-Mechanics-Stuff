%% This code calculates topocentric right ascension and declination for ISS, given its geocentric Latitude, Longitude & Altitude at time instant
clear
close all
clc

%% Dimensions in km
Re = 6378;
f = 0.003353;
theta_G = 126.7; % [deg]
r = [-5368; -1784; 3691];

%% Topocentric Origin (Observer Location)
H = 0;
Lat = 20; % [deg]
Lon = 60; % [deg]

theta = theta_G + Lon;
R = [(Re/sqrt(1-(2*f-f^2)*sind(Lat)^2)+H)*cosd(Lat)*cosd(theta); (Re/sqrt(1-(2*f-f^2)*sind(Lat)^2)+H)*cosd(Lat)*sind(theta); (Re*(1-f)^2/sqrt(1-(2*f-f^2)*sind(Lat)^2)+H)*sind(Lat)];
rho_vec = r - R;
rho = norm(rho_vec);
delta = asind(rho_vec(3)/rho);
alpha = asind(rho_vec(2)/rho/cosd(delta));

