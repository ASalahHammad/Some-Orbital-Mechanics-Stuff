clear
clc
close all

G = 6.6743e-11;

% Sun-Earth
tf = 2*365*24*3600; % [s]
dt = 100;
t = 0:dt:tf;
m1 = 1.989e30;
m2 = 5.974e24;
r1_0 = [0; 0; 0];
v1_0 = [0e3; 0e3; 0e3];
r2_0 = [149.6e9; 0; 0];
v2_0 = [0; 29.78e3; 0];
obj1.filename = "sun.png";
obj1.radius = 696e7;
obj2.filename = "earth.png";
obj2.radius = .4e10;

% Earth-Moon
% tf = 2*30*24*3600; % [s]
% dt = 10;
% t = 0:dt:tf;
% m1 = 5.974e24;
% m2 = 73.48e21;
% r1_0 = [0; 0; 0];
% v1_0 = [0; 0; 0];
% r2_0 = [384.4e6; 0; 0];
% v2_0 = [0; 1.022e3; 0];
% obj1.filename = "earth.png";
% obj1.radius = 3e7;
% obj2.filename = "moon.png";
% obj2.radius = 2e7;

%% Solve
[t, states] = RK4(@(t, states_0) fdot(t, states_0, m1, m2, G), t, [r1_0; r2_0; v1_0; v2_0]);

X1 = states(:, 1).';
Y1 = states(:, 2).';
Z1 = states(:, 3).';
X2 = states(:, 4).';
Y2 = states(:, 5).';
Z2 = states(:, 6).';

R1 = [X1; Y1; Z1];
R2 = [X2; Y2; Z2];
RG = (m1*R1 + m2*R2) / (m1 + m2);
V1 = states(:, 7:9).';
V2 = states(:, 10:12).';
VG = (m1*V1 + m2*V2) / (m1 + m2);
r = sqrt((X2-X1).^2 + (Y2-Y1).^2 + (Z2-Z1).^2);
A1 = G*m2*(R2-R1)./r.^3;
A2 = -G*m1*(R2-R1)./r.^3;
AG = (m1*A1 + m2*A2) / (m1 + m2);

%% Animate & Plot
fig1 = figure('Position', [1, 1, 1366, 728]);
ax1 = gca;
% xlim([0, 6*10^6]);
% ylim([0, 14e6]);
% zlim([0, 7e6]);
simulate(t, 20000, [R1; R2; RG], [V1; V2; VG], [A1; A2; AG], ["r", "g", "b--"], ["body1", "body2", "CG"], ax1, [obj1, obj2]);

%%
% fig2 = figure; ax2 = gca;
% % xlim([-2e6, 4*10^6]);
% % ylim([-2e6, 5e6]);
% % zlim([-2e6, 2e6]);
% simulate(t, 2000, [R1; R2; RG]-[RG;RG;RG], [V1; V2; VG], [A1; A2; AG], ["r", "g", "b--"], ["body1", "body2", "CG"], ax2);

%%
% fig3 = figure; ax3 = gca;
% % xlim([-2e6, 4*10^6]);
% % ylim([-2e6, 5e6]);
% % zlim([-2e6, 2e6]);
% simulate(t, 2000, [R1; R2; RG]-[R1;R1;R1], [V1; V2; VG], [A1; A2; AG], ["r", "g", "b--"], ["body1", "body2", "CG"], ax3);

%%
% fig4 = figure; ax3 = gca;
% % xlim([-5e6, 4*10^6]);
% % ylim([-2e6, 5e6]);
% % zlim([-2e6, 2e6]);
% simulate(t, 2000, [R1; R2; RG]-[R2;R2;R2], [V1; V2; VG], [A1; A2; AG], ["r", "g", "b--"], ["body1", "body2", "CG"], ax3);

%%
V1_mag = vecnorm(V1);
A1_mag = vecnorm(A1);

figure;
hold on; grid on; 
xlabel("time"); ylabel("$||V||$", 'Interpreter', 'latex');
plot(t, V1_mag);
title("Velocity of Body 1");

figure;
hold on; grid on;
xlabel("time"); ylabel("$||A||$", 'Interpreter', 'latex');
plot(t, A1_mag);
title("Acceleration of Body 1");


%% Plot V, A, Sun-Earth-Sat problem, plot velocity distribution over orbit
