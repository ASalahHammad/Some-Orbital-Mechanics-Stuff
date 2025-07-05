clear
close all
clc

r0 = 10000; % [km]
v0 = 10; % [km/s]
theta0 = 30*pi/180; % [rad]
mu = 398600; % [km^3/s^2]
e = roots([1, 2*cos(theta0)-cos(theta0)*r0*v0^2/mu, 1-r0*v0^2/mu]);
e = e(e>=0);
h = sqrt(r0*mu*(1+e*cos(theta0)));
a = h^2/mu/(1-e^2);
alpha = 1/a;
vr0 = mu/h*e*sin(theta0);

%% Find Universal Anomaly
dt = 3600; % [s] fix this
S = @(z) [((sinh(sqrt(-z(z<0)))-sqrt(-z(z<0)))./sqrt(-z(z<0)).^3), 1/6*z(z==0)^0, ((sqrt(z(z>0))-sin(sqrt(z(z>0))))./sqrt(z(z>0)).^3)];
C = @(z) [(cosh(sqrt(-z(z<0)))-1)./-z(z<0), 1/2*z(z==0)^0, (1-cos(sqrt(z(z>0))))./z(z>0)];
figure;
subplot(1,3,1);
hold on; grid on;
z = -50:0;
plot(z, S(z)); hold on; grid on; plot(z, C(z));
legend("S(z)", "C(z)");
subplot(1,3,2);
hold on; grid on;
z = 0:30;
plot(z, S(z)); hold on; grid on; plot(z, C(z));
legend("S(z)", "C(z)");
subplot(1,3,3);
hold on; grid on;
z = 0:500;
plot(z, S(z)); hold on; grid on; plot(z, C(z));
legend("S(z)", "C(z)");
ylim([0,0.04]);

S_d = @(z) 1/(2*z) * (C(z)-3*S(z));
C_d = @(z) 1/(2*z) * (1-z*S(z)-2*C(z));



f = @(Kai) r0*vr0/sqrt(mu)*Kai^2*C(alpha*Kai^2) + (1-alpha*r0)*Kai^3*S(alpha*Kai^2) + r0*Kai - sqrt(mu)*dt;
f_d = @(Kai) 2*r0*vr0/sqrt(mu)*Kai*C(alpha*Kai^2) + r0*vr0/sqrt(mu)*Kai^2*C_d(alpha*Kai^2)*(2*alpha*Kai) + 3*(1-alpha*r0)*Kai^2*S(alpha*Kai^2) + (1-r0*alpha)*Kai^3*S_d(alpha*Kai^2)*(2*alpha*Kai) + r0;
Kai = newton_raphson(f, f_d, sqrt(mu)*abs(alpha)*dt, 1e-6)
F0 = 2*atanh(sqrt((e-1)/(e+1))*tan(theta0/2))
F = F0 + Kai/sqrt(-a)
theta = 2*atan(sqrt((e+1)/(e-1))*tanh(F/2))
disp(theta*180/pi)

