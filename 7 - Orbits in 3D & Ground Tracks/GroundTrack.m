clear
close all
clc

% X means Inertial frame
% x_dash means ECEF frame
% x_hat means in the perifocal frame

%% Video File
videofile = VideoWriter("Ground Track.avi");
open(videofile);

%% Problem Data
mu = 398600; % [km^3/s^2]
J2 = 1.08263e-03;
R_earth = 6378; % [km]s
rp  = 6700; % [km]
ra = 10000; % [km]
theta0 = 230 *pi/180; % [rad]
OMEGA0 = 270 *pi/180; % [rad]
i0 = 60 *pi/180; % [rad]
omega0 = 45 *pi/180; % [rad]
delta_t = 1*60; % [s]

a = (rp+ra)/2;
e = (ra-rp)/(ra+rp);
h = sqrt(mu*rp*(1+e)); % [km^2/s]
T = 2*pi/sqrt(mu)*a^(3/2); % [s]
OMEGA_d = -1.5*sqrt(mu)*J2*R_earth^2/(1-e^2)^2/a^(7/2) * cos(i0); % [rad/s]
% omega_d = -1.5*sqrt(mu)*J2*R_earth^2/(1-e^2)^2/a^(7/2) * (5/2*sind(i0)^2-2)
omega_d = OMEGA_d * (5/2*sin(i0)^2-2)/cos(i0); % [rad/s]

E0 = 2*atan(tan(theta0/2) * sqrt((1-e)/(1+e))); % [rad]
M0 = E0 - e*sin(E0); % [rad]
t0 = M0/2/pi*T;

N_points = 412;
delta = nan(1,N_points);
alpha = nan(1,N_points);
r_X = nan(3, N_points);
r_xdash = nan(3, N_points);

for n_time=1:N_points
    dt = (n_time-1)*delta_t;
    t = t0 + dt;
    M = 2*pi*t/T;
    E = newton_raphson(@(E) E-e*sin(E)-M, @(E) 1-e*cos(E), E0); % this might need fixing for E to become +ve
    theta = 2*atan(tan(E/2) * sqrt((1+e)/(1-e))); % [rad]
    OMEGA = OMEGA0 + OMEGA_d*dt;
    omega = omega0 + omega_d*dt;
    i = i0;

    r_x_hat = h^2/mu/(1+e*cos(theta)) * [cos(theta); sin(theta); 0]; % in orbit axes
    v_x_hat = mu/h * [-sin(theta); e+cos(theta); 0];
    Q_X_xhat = [cos(omega),sin(omega),0; -sin(omega),cos(omega),0; 0,0,1] * [1,0,0; 0,cos(i),sin(i); 0,-sin(i),cos(i)] * [cos(OMEGA),sin(OMEGA),0; -sin(OMEGA),cos(OMEGA),0; 0,0,1];
    r_X(:,n_time) = transpose(Q_X_xhat) * r_x_hat; % in inertial axes fixed with respect to stars
    v_X = transpose(Q_X_xhat) * v_x_hat;

    OMEGA_EARTH = 2*pi/24/3600; % [rad/s] anglar velocity of earth about inertial Z axis
    THETA = OMEGA_EARTH*dt;
    r_xdash(:, n_time) = [cos(THETA),sin(THETA),0; -sin(THETA),cos(THETA),0; 0,0,1] * r_X(:,n_time); % in earth fixed axes

    r_xdash_mag = norm(r_xdash(:,n_time));
    l = r_xdash(1,n_time)/r_xdash_mag;
    m = r_xdash(2,n_time)/r_xdash_mag;
    ind = r_xdash(3,n_time)/r_xdash_mag;

    delta(n_time) = asin(ind); % [rad] declination
    if(m>0), alpha(n_time) = acos(l/cos(delta(n_time))); else, alpha(n_time) = 2*pi-acos(l/cos(delta(n_time))); end  % [rad] right ascension
end

alpha(alpha>pi) = alpha(alpha>pi)-2*pi; % this is to adjust the vector for plotting

%% Plot
% fig = figure('Position', [1, 1, 1368, 728]);
% plot(ax1, alpha*180/pi, delta*180/pi, 'rx', 'LineWidth', 2);

fig = figure('Position', [1, 1, 1368, 728]);
ax1 = subplot(2,15,[1:10,16:25]);
IMAGE = imread("earth.png");
X = -180:180;
Y = 90:-1:-90;
p1 = plot(ax1, nan, nan);
imagesc(X, Y, IMAGE);
axis xy;
hold on; grid on;
title("Ground Track of Satellite", 'Interpreter', 'latex', 'FontSize', 15);
xlabel("Right Ascension $\alpha$ [deg]", 'Interpreter', 'latex', 'FontSize', 15);
ylabel("Declination $\delta$ [deg]", 'Interpreter', 'latex', 'FontSize', 15);

ax2 = subplot(2,15,12:15);
hold on; grid on;
p2 = plot3(ax2, nan, nan, nan, nan, nan, nan);
p3 = plot3(ax2, nan, nan, nan, nan, nan, nan);
[X_sp1, Y_sp1, Z_sp1] = sphere;
X_sp1 = X_sp1*R_earth; Y_sp1 = Y_sp1*R_earth; Z_sp1 = -Z_sp1*R_earth;
surf(ax2, X_sp1, Y_sp1, Z_sp1, 'CData', IMAGE, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
title("Satellite Orbit in ECEF axes", 'Interpreter', 'latex', 'FontSize', 15);
xlabel("$x$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
ylabel("$y$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
zlabel("$z$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
view([150,10])
xlim([-6000, 6000]);
ylim([-1e4, 1e4]);
zlim([-1e4, 1e4]);
axis equal;

ax3 = subplot(2,15,27:30);
hold on; grid on;
[X_sp2, Y_sp2, Z_sp2] = sphere;
X_sp2 = X_sp2*R_earth; Y_sp2 = Y_sp2*R_earth; Z_sp2 = -Z_sp2*R_earth;
p4 = plot3(ax2, nan, nan, nan, nan, nan, nan);
p5 = plot3(ax2, nan, nan, nan, nan, nan, nan);
p6 = surf(ax3, X_sp2, Y_sp2, Z_sp2, 'CData', IMAGE, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
title("Satellite Orbit in Inertial axes", 'Interpreter', 'latex', 'FontSize', 15);
xlabel("$x$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
ylabel("$y$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
zlabel("$z$ [km]", 'Interpreter', 'latex', 'FontSize', 12);
view([150,10])
xlim([-6000, 6000]);
ylim([-1e4, 1e4]);
zlim([-1e4, 1e4]);
axis equal;

for ind = 1:N_points
    % delete(p1);
    delete(p2);
    delete(p3);
    delete(p4);
    delete(p5);
    delete(p6);
    p1 = plot(ax1, alpha(1:ind)*180/pi, delta(1:ind)*180/pi, 'rx', 'LineWidth', 2);
    p2 = plot3(ax2, r_xdash(1,1:ind), r_xdash(2,1:ind), r_xdash(3,1:ind), 'r', 'LineWidth', 2);
    p3 = plot3(ax2, r_xdash(1,ind), r_xdash(2,ind), r_xdash(3,ind), 'ro', 'LineWidth', 4);
    X_sp2_temp = cos(OMEGA_EARTH*delta_t)*X_sp2 - sin(OMEGA_EARTH*delta_t)*Y_sp2;
    Y_sp2 = sin(OMEGA_EARTH*delta_t)*X_sp2 + cos(OMEGA_EARTH*delta_t)*Y_sp2;
    X_sp2 = X_sp2_temp;
    clear X_sp2_temp;
    p4 = plot3(ax3, r_X(1,1:ind), r_X(2,1:ind), r_X(3,1:ind), 'r', 'LineWidth', 2);
    p5 = plot3(ax3, r_X(1,ind), r_X(2,ind), r_X(3,ind), 'ro', 'LineWidth', 4);
    p6 = surf(ax3, X_sp2, Y_sp2, Z_sp2, 'CData', IMAGE, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    pause(0.01);
    frame = getframe(gcf);
    writeVideo(videofile,frame);
end

close(videofile);

saveas(fig, "Ground Track", 'svg');

%% NACA Orbit Parameters, 3d simulation along with 2d track simulation
