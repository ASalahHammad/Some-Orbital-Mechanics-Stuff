clear
close all
clc

m1 = 5.974e24; % earth
m2 = 0; % sat
G = 6.6743e-11;
mu = G*(m1+m2);

r0_vec = [7000e3; 0];
r0 = norm(r0_vec);
theta0 = 0;
e = 0.2;

h = sqrt(r0*mu*(1+e*cos(theta0)));

v0 = h/r0;
v0_vec = [h/(r0_vec(1)*tan(theta0+pi/2) - r0_vec(2)); h*tan(theta0+pi/2)/(r0_vec(1)*tan(theta0+pi/2) - r0_vec(2))];

f_func = @(dt) 1 - mu/2/r0^3*dt.^2 + mu/2*dot(r0_vec,v0_vec)/r0^5*dt.^3 + mu/24*(-2*mu/r0^6+3*v0^2/r0^5-15*dot(r0_vec,v0_vec)^2/r0^7)*dt.^4;
g_func = @(dt) dt - 1/6*mu/r0^3*dt.^3 + mu/4*dot(r0_vec,v0_vec)/r0^5*dt.^4;

dt = 0:960;
f = f_func(dt);
g = g_func(dt);

r_vec = f.*r0_vec + g.*v0_vec;

figure;
grid minor; hold on;
plot(dt, vecnorm(r_vec, 2, 1)/1000, 'k', 'LineWidth', 2);
title("Comparison between exact solution and the $f$ \& $g$ series", 'interpreter', 'latex');
xlabel("time [s]", 'Interpreter', 'latex'); ylabel("r [km]", 'Interpreter', 'latex');
legend("f and g series");
