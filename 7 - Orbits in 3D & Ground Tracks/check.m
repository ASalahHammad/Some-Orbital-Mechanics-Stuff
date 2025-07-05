clc
close all
clear

mu = 398600; % [km^3/s^2]
R_earth = 6378; % [km]
e = 0.001;
J2 = 1.08263e-03;
a = R_earth+(300:200:1100).'; % [km]

i_vec = 0:5:100;

OMEGA_d_vec = -3/2*sqrt(mu)*J2*R_earth^2./((1-e^2)^2*a.^(7/2)).*cosd(i_vec)           *3600*24*180/pi;
omega_d_vec = -3/2*sqrt(mu)*J2*R_earth^2./((1-e^2)^2*a.^(7/2)).*(5/2*sind(i_vec).^2-2)*3600*24*180/pi;

figure('Position', [1, 1, 1368, 728]);
subplot(1,2,1);
plot(i_vec,OMEGA_d_vec,'LineWidth',2);
grid on; xlabel('$i$, degrees','Interpreter','latex','FontSize',15); ylabel('$\dot{\Omega}$, degrees per day','Interpreter','latex','FontSize',15);
legend(num2str(a),'Location','best');
subplot(1,2,2);
plot(i_vec,omega_d_vec,'LineWidth',2);
grid on; xlabel('$i$, degrees','Interpreter','latex','FontSize',15); ylabel('$\dot{\omega}$, degrees per day','Interpreter','latex','FontSize',15);
legend(num2str(a),'Location','best');
