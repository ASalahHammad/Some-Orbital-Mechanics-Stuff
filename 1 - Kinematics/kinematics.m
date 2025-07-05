close all
clear
clc

tf = 1.5;
dt = 0.02;
t_vec = 0:dt:tf;

a = [16*ones(size(t_vec)); 30*t_vec; 18/5*t_vec.^2+4];
r0 = [6; 4; 1];
v0 = [7; 0; 0];

% [t, states] = ode45(@fdot, t_vec, [r0, v0]);
[t, states] = RK4(@fdot, t_vec, [r0; v0]);
r = states(:, 1:3).'; % I need r & v as 3xn rows
v = states(:, 4:6).';

simulate(t_vec, r, v, a);
