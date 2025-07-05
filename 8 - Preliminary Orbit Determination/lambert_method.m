%% Algorithm 5.2
clear
close all
clc

%% All dimensions in km
mu = 398600;
Re = 6378;

r1_vec = [5000; 10000; 2100];
r2_vec = [-14600; 2500; 7000];
dt = 3600; % [s]

r1 = norm(r1_vec);
r2 = norm(r2_vec);

%% Assume a prograde trajectory
h_hat = cross(r1_vec, r2_vec) / r1 / r2;
if(h_hat(3)>=0)
    dtheta = acosd(dot(r1_vec,r2_vec)/r1/r2); % [rad]
else
    dtheta = 360 - acosd(dot(r1_vec,r2_vec)/r1/r2); % [rad]
end
A = sind(dtheta) * sqrt(r1*r2/(1-cosd(dtheta)));
%% Solve for z
S = @(z) [(sinh(sqrt(-z(z<0)))-sqrt(-z(z<0)))./sqrt(-z(z<0)).^3, 1/6*z(z==0)^0, (sqrt(z(z>0))-sin(sqrt(z(z>0))))./sqrt(z(z>0).^3)];
C = @(z) [(cosh(sqrt(-z(z<0)))-1)./-z(z<0), 1/2*z(z==0)^0, (1-cos(sqrt(z(z>0))))./z(z>0)];
S_d = @(z) 1/(2*z) * (C(z)-3*S(z));
C_d = @(z) 1/(2*z) * (1-z*S(z)-2*C(z));
y = @(z) r1 + r2 + A*(z*S(z)-1)/sqrt(C(z));
y_d = @(z) A/4*sqrt(C(z));
F = @(z) (y(z)/C(z))^(3/2)*S(z) + A*sqrt(y(z)) - sqrt(mu)*dt;

F_d = @(z) feval(@(x) x{1+(z~=0)}, {sqrt(2)/40*y(0)^(3/2) + A/8*(sqrt(y(0) + A*sqrt(1/2/y(0)))),     (y(z)/C(z))^(3/2)*(1/(2*z)*(C(z)-1.5*S(z)/C(z))+0.75*S(z)^2/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z))+A*sqrt(C(z)/y(z)))}); %#ok

Z = newton_raphson(F, F_d, 1.5);

%% Get lagrange coefficients
f = 1-y(Z)/r1;
g = A*sqrt(y(Z)/mu);
f_d = sqrt(mu)/r1/r2*sqrt(y(Z)/C(Z))*(Z*S(Z)-1);
g_d = 1-y(Z)/r2;

v1_vec = 1/g*(r2_vec-f*r1_vec);
v2_vec = 1/g*(g_d*r2_vec-r1_vec);
v1 = norm(v1_vec);
v2 = norm(v2_vec);

%% You can now use r1 & v1 or r2 & v2 to proceed with algorithm 4.2
[a,e,i,OMEGA,omega,theta] = orbital_elements(r2_vec,v2_vec)
h = norm(cross(r2_vec, v2_vec));
rp = h^2/mu/(1+e);
T = 2*pi/sqrt(mu)*a^(3/2) / 3600; % hr
% perigee altitude
perigee_alt = rp - Re;
% time since perigee
E1 = 2*atan(sqrt((1-e)/(1+e))*tand(theta/2)); % rad
Me1 = E1 - e*sin(E1);
t1 = h^3/mu^2/(1-e^2)^(3/2)*Me1;
