%% Algorithms 5.5 & 5.6
clear
clc
close all
format long

%% All dimensions in km, angles in deg, and angular rates in rad/s
mu = 398600;
f = 1/298.257223563; % Flattening of the Ellipsoid
Re = 6378; % [Km] Earth's Semi Major Axis, sometimes called a

%% Observer Location
phi = 40; % [deg]
H = 1; % [Km]

%% reading 1
t1 = 0;
alpha1 = 43.537; % [deg]
delta1 = -8.7833; % [deg]
theta1 = 44.506; % [deg]
R1_vec = (Re/sqrt(1-(2*f-f^2)*sind(phi)^2)+H) * [cosd(phi)*cosd(theta1); cosd(phi)*sind(theta1); (1-f)^2*sind(phi)];
rho_hat_vec1 = [cosd(delta1)*cosd(alpha1); cosd(delta1)*sind(alpha1); sind(delta1)];
%% reading 2
t2 = 118.1;
alpha2 = 54.420; % [deg]
delta2 = -12.074; % [deg]
theta2 = 45.000; % [deg]
R2_vec = (Re/sqrt(1-(2*f-f^2)*sind(phi)^2)+H) * [cosd(phi)*cosd(theta2); cosd(phi)*sind(theta2); (1-f)^2*sind(phi)];
rho_hat_vec2 = [cosd(delta2)*cosd(alpha2); cosd(delta2)*sind(alpha2); sind(delta2)];
%% reading 3
t3 = 237.58;
alpha3 = 64.318; % [deg]
delta3 = -15.105; % [deg]
theta3 = 45.499; % [deg]
R3_vec = (Re/sqrt(1-(2*f-f^2)*sind(phi)^2)+H) * [cosd(phi)*cosd(theta3); cosd(phi)*sind(theta3); (1-f)^2*sind(phi)];
rho_hat_vec3 = [cosd(delta3)*cosd(alpha3); cosd(delta3)*sind(alpha3); sind(delta3)];

%% Algorithm 5.5
tau1 = t1 - t2;
tau3 = t3 - t2;
tau = t3 - t1;
p1_vec = cross(rho_hat_vec2,rho_hat_vec3);
p2_vec = cross(rho_hat_vec1,rho_hat_vec3);
p3_vec = cross(rho_hat_vec1,rho_hat_vec2);
D0 = dot(rho_hat_vec1,p1_vec);
D11 = dot(R1_vec,p1_vec);
D21 = dot(R2_vec,p1_vec);
D31 = dot(R3_vec,p1_vec);
D12 = dot(R1_vec,p2_vec);
D22 = dot(R2_vec,p2_vec);
D32 = dot(R3_vec,p2_vec);
D13 = dot(R1_vec,p3_vec);
D23 = dot(R2_vec,p3_vec);
D33 = dot(R3_vec,p3_vec);
A = 1/D0*(-D12*tau3/tau+D22+D32*tau1/tau);
B = 1/6/D0 * (D12*(tau3^2-tau^2)*tau3/tau + D32*(tau^2-tau1^2)*tau1/tau);
E = dot(R2_vec,rho_hat_vec2);
a = -(A^2+2*A*E+norm(R2_vec)^2);
b = -2*mu*B*(A+E);
c = -mu^2*B^2;
r2 = newton_raphson(@(x) x^8+a*x^6+b*x^3+c, @(x) 8*x^7+6*a*x^5+3*b*x^2, 9000); % this result should be dealt with very carefully
rho1 = 1/D0*((6*(D31*tau1/tau3+D21*tau/tau3)*r2^3+mu*D31*(tau^2-tau1^2)*tau1/tau3) / (6*r2^3+mu*(tau^2-tau3^2)) - D11);
rho2 = A + mu*B/r2^3;
rho3 = 1/D0*((6*(D13*tau3/tau1-D23*tau/tau1)*r2^3+mu*D13*(tau^2-tau3^2)*tau3/tau1) / (6*r2^3+mu*(tau^2-tau1^2)) - D33);
r1_vec = R1_vec + rho1*rho_hat_vec1;
r2_vec = R2_vec + rho2*rho_hat_vec2;
r3_vec = R3_vec + rho3*rho_hat_vec3;
f1 = 1 - 0.5*mu/r2^3*tau1^2;
f3 = 1 - 0.5*mu/r2^3*tau3^2;
g1 = tau1 - 1/6*mu/r2^3*tau1^3;
g3 = tau3 - 1/6*mu/r2^3*tau3^3;
v2_vec = 1/(f1*g3-f3*g1)*(-f3*r1_vec+f1*r3_vec);

r_vec = r2_vec;
v_vec = v2_vec;
r = norm(r2_vec);
v = norm(v2_vec);
vr = dot(v2_vec,r_vec/r);

old_rho = [rho1,rho2,rho3];
epsilon = 1;
%% Algorithm 5.6 Iterative improvement of the orbit
% while(epsilon>1e-5)
%     alpha = 2/r-v^2/mu;
%     % Algorithm 3.3 to find universal variable Kai
%     S = @(z) [((sinh(sqrt(-z(z<0)))-sqrt(-z(z<0)))./sqrt(-z(z<0)).^3), 1/6*z(z==0)^0, ((sqrt(z(z>0))-sin(sqrt(z(z>0))))./sqrt(z(z>0)).^3)];
%     C = @(z) [(cosh(sqrt(-z(z<0)))-1)./-z(z<0), 1/2*z(z==0)^0, (1-cos(sqrt(z(z>0))))./z(z>0)];
%     S_d = @(z) 1/(2*z) * (C(z)-3*S(z));
%     C_d = @(z) 1/(2*z) * (1-z*S(z)-2*C(z));
%     f = @(Kai,dt) r*vr/sqrt(mu)*Kai^2*C(alpha*Kai^2) + (1-alpha*r)*Kai^3*S(alpha*Kai^2) + r*Kai - sqrt(mu)*dt;
%     f_d = @(Kai) 2*r*vr/sqrt(mu)*Kai*C(alpha*Kai^2) + r*vr/sqrt(mu)*Kai^2*C_d(alpha*Kai^2)*(2*alpha*Kai) + 3*(1-alpha*r)*Kai^2*S(alpha*Kai^2) + (1-r*alpha)*Kai^3*S_d(alpha*Kai^2)*(2*alpha*Kai) + r;
%     Kai1 = newton_raphson(@(Kai) f(Kai,tau1), f_d, sqrt(mu)*abs(alpha)*tau1);
%     f1 = (1-Kai1^2/r*C(alpha*Kai1^2) + f1)/2;
%     g1 = (tau1 - 1/sqrt(mu)*Kai1^3*S(alpha*Kai1^2) + g1)/2;
%     Kai3 = newton_raphson(@(Kai) f(Kai,tau3), f_d, sqrt(mu)*abs(alpha)*tau3);
%     f3 = (1-Kai3^2/r*C(alpha*Kai3^2) + f3)/2;
%     g3 = (tau3 - 1/sqrt(mu)*Kai3^3*S(alpha*Kai3^2) + g3)/2;
%     c1 = g3/(f1*g3-f3*g1);
%     c3 = -g1/(f1*g3-f3*g1);
%     rho1 = 1/D0*((6*(D31*tau1/tau3+D21*tau/tau3)*r2^3+mu*D31*(tau^2-tau1^2)*tau1/tau3) / (6*r2^3+mu*(tau^2-tau3^2)) - D11);
%     rho2 = A + mu*B/r2^3;
%     rho3 = 1/D0*((6*(D13*tau3/tau1-D23*tau/tau1)*r2^3+mu*D13*(tau^2-tau3^2)*tau3/tau1) / (6*r2^3+mu*(tau^2-tau1^2)) - D33);
%     r1_vec = R1_vec + rho1*rho_hat_vec1;
%     r2_vec = R2_vec + rho2*rho_hat_vec2;
%     r3_vec = R3_vec + rho3*rho_hat_vec3;
%     v2_vec = 1/(f1*g3-f3*g1)*(-f3*r1_vec+f1*r3_vec);
%     r_vec = r2_vec;
%     v_vec = v2_vec;
%     r = norm(r2_vec);
%     v = norm(v2_vec);
%     vr = dot(v2_vec,r_vec/r);
%     % check
%     epsilon = max([rho1,rho2,rho3]-old_rho);
%     old_rho = [rho1,rho2,rho3];
% end

%% Algorithm 4.2 to get orbital elements
[a,e,i,OMEGA,omega,theta] = orbital_elements(r_vec,v_vec)
