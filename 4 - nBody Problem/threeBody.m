clear
clc
close all

G = 6.6743e-11;

%% Book's Problem
tfinal = 67000;
dt = 5;
t = 0:dt:tfinal;
m1 = 1.e29;
m2 = 1.e29;
m3 = 1.e29;
r1_0 = [0; 0; 0];
v1_0 = [0e3; 0e3; 0e3];
r2_0 = [3.e8; 0; 0];
v2_0 = [250e3; 250e3; 0];
r3_0 = [6.e8; 0; 0];
v3_0 = [0; 0e3; 0];
obj1.filename = "./Images/sun.png";
obj1.radius = 696e5;
obj2.filename = "./Images/earth.png";
obj2.radius = .4e8;
obj3.filename = "./Images/venus.jpg";
obj3.radius = .4e8;

% Sun-Earth
% tfinal = 2*365*24*3600; % [s]
% dt = 100;
% t = 0:dt:tfinal;
% m1 = 1.989e30;
% m2 = 5.974e24;
% r1_0 = [0; 0; 0];
% v1_0 = [0e3; 0e3; 0e3];
% r2_0 = [149.6e9; 0; 0];
% v2_0 = [0; 29.78e3; 0];
% obj1.filename = "sun.png";
% obj1.radius = 696e7;
% obj2.filename = "earth.png";
% obj2.radius = .4e10;

% Earth-Moon
% tfinal = 2*30*24*3600; % [s]
% dt = 10;
% t = 0:dt:tfinal;
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
[t, states] = RK4(@(t, states_0) fdot(t, states_0, [m1, m2, m3]), t, [r1_0; r2_0; r3_0; v1_0; v2_0; v3_0]);

R1 = states(:, 1:3).';
R2 = states(:, 4:6).';
R3 = states(:, 7:9).';
RG = (m1*R1 + m2*R2 + m3*R3) / (m1 + m2 + m3);
V1 = states(:, 10:12).';
V2 = states(:, 13:15).';
V3 = states(:, 16:18).';
VG = (m1*V1 + m2*V2 + m3*V3) / (m1 + m2 + m3);

r12 = norm(R2-R1);
r13 = norm(R3-R1);
r23 = norm(R3-R2);

A1 = G*m2*(R2-R1)./r12.^3 + G*m3*(R3-R1)./r13.^3;
A2 = G*m1*(R1-R2)./r12.^3 + G*m3*(R3-R2)./r23.^3;
A3 = G*m1*(R1-R3)./r13.^3 + G*m2*(R2-R3)./r23.^3;
AG = (m1*A1 + m2*A2+ m3*A3) / (m1 + m2+ m3);

%% Animate & Plot
fig1 = figure('Position', [1, 1, 1366, 728]);
ax1 = gca;
% axis equal
hold on;
grid on;
% xlim([0, 6*10^6]);
% ylim([0, 14e6]);
% zlim([0, 7e6]);
% view(2);
xlabel("x [Km]"); ylabel("y [Km]"); zlabel("z [Km]");
simulate(t, 200, [R1-RG(:,1); R2-RG(:,1); R3-RG(:,1); RG-RG(:,1)], [V1-VG(:,1); V2; V3; VG], [A1; A2; A3; AG], ["r", "g", "b", "k--"], ["Sun", "Earth", "venus", "CG"], ax1, [obj1, obj2, obj3]);

% fig2 = figure('Position', [1, 1, 1366, 728]);
% ax2 = gca;
% axis equal
% hold on;
% grid on;
% % pbaspect([1,1,1]);
% xlim([-1.2e9, 1.5e9]);
% ylim([-5.1e8, 8e8]);
% % zlim([0, 7e6]);
% view(2);
% xlabel("x [Km]"); ylabel("y [Km]"); zlabel("z [Km]");
% simulate(t, 200, [R1-RG; R2-RG; R3-RG; RG-RG], [V1-VG; V2-VG; V3-VG; VG-VG], [A1; A2; A3; AG], ["r", "g", "b", "k--"], ["Sun", "Earth", "venus", "CG"], ax2, [obj1, obj2, obj3]);
