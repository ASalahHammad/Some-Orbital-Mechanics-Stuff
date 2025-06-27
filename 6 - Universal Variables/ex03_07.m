clear
close all
clc

mu = 398600; % [km^3/s^2]

r0_vec = [7000; -12124]; % [km]
v0_vec = [2.6679; 4.6210]; % [km/s]

r0 = norm(r0_vec)
v0 = norm(v0_vec)

vr0 = dot(r0_vec,v0_vec) / r0
alpha = 2/r0 - v0^2/mu

%% Algorithm 3.3
dt = 3600; % [s] fix this
S = @(z) [(sinh(sqrt(-z(z<0)))-sqrt(-z(z<0)))./sqrt(-z(z<0)).^3, 1/6*z(z==0)^0, (sqrt(z(z>0))-sin(sqrt(z(z>0))))./sqrt(z(z>0).^3)];
C = @(z) [(cosh(sqrt(-z(z<0)))-1)./-z(z<0), 1/2*z(z==0)^0, (1-cos(sqrt(z(z>0))))./z(z>0)];
S_d = @(z) 1/(2*z) * (C(z)-3*S(z));
C_d = @(z) 1/(2*z) * (1-z*S(z)-2*C(z));
func = @(Kai) r0*vr0/sqrt(mu)*Kai^2*C(alpha*Kai^2) + (1-alpha*r0)*Kai^3*S(alpha*Kai^2) + r0*Kai - sqrt(mu)*dt;
func_d = @(Kai) 2*r0*vr0/sqrt(mu)*Kai*C(alpha*Kai^2) + r0*vr0/sqrt(mu)*Kai^2*C_d(alpha*Kai^2)*(2*alpha*Kai) + 3*(1-alpha*r0)*Kai^2*S(alpha*Kai^2) + (1-r0*alpha)*Kai^3*S_d(alpha*Kai^2)*(2*alpha*Kai) + r0;
Kai = newton_raphson(func, func_d, sqrt(mu)*abs(alpha)*dt, 1e-6)
% Kai = 253.53 % this is weird and needs fixing

%% lagrange coefficients, r, v
f = 1-Kai^2/r0*C(alpha*Kai^2)
g = dt - 1/sqrt(mu)*Kai^3*S(alpha*Kai^2)

r_vec = f*r0_vec + g*v0_vec
r = norm(r_vec)

fd = sqrt(mu)/r/r0 * (alpha*Kai^3*S(alpha*Kai^2) - Kai)
gd = 1 - Kai^2/r*C(alpha*Kai^2)

v_vec = fd*r0_vec + gd*v0_vec
v = norm(v_vec)
