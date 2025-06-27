%% =============== Algorithm 2.3 ===================

clear
close all
clc

global r_func f_func g_func f_d_func g_d_func x0 y0 r0_vec v0_vec r_vec_traj r_max_pos

m1 = 5.974e24; % earth
m2 = 0; % sat
G = 6.6743e-11;
mu = G*(m1+m2);

x0 = 8182.4e3;
y0 = -6865.9e3;
x0_d = 0.47572e3;
y0_d = 8.8116e3;
r0_vec = [x0; y0];
v0_vec = [x0_d; y0_d];

r0 = sqrt(x0^2 + y0^2);
v0 = sqrt(x0_d^2 + y0_d^2);

vr0 = (x0*x0_d + y0*y0_d) / r0; % radial component of v0, i.e. projection of v0 on r0
h = r0*sqrt(v0^2 - vr0^2);

r_func = @(dtheta) h^2/mu ./ (1 + (h^2/mu/r0-1)*cos(dtheta) - h*vr0/mu*sin(dtheta));

% lagrange coefficients
f_func = @(dtheta) 1 - mu*r_func(dtheta)/h^2.*(1-cos(dtheta));
g_func = @(dtheta) r_func(dtheta)*r0/h.*sin(dtheta);
f_d_func = @(dtheta) mu/h*(1-cos(dtheta))./sin(dtheta) .* (mu/h^2*(1-cos(dtheta)) - 1/r0 - 1./r_func(dtheta));
g_d_func = @(dtheta) 1-mu*r0/h^2*(1-cos(dtheta));

%% Total Solution
dtheta_vec = linspace(-50, 120, 50) * pi/180;
r = r_func(dtheta_vec);
f = f_func(dtheta_vec);
g = g_func(dtheta_vec);
f_d = f_d_func(dtheta_vec);
g_d = g_d_func(dtheta_vec);

r_vec_traj = f.*r0_vec + g.*v0_vec;
v_vec_traj = f_d.*r0_vec + g_d.*v0_vec;

e = sqrt((h^2/r0/mu-1)^2 + (vr0*h/mu)^2);
theta0 = atan2(vr0*h*r0, h^2-mu*r0);

% Plot
fig = figure('Position', [1, 1, 1366, 650]);
temp.Value = 120;
sliderCallback(temp);
clear temp;
slider = uicontrol("Style","slider", "Min",-50, "Max",120, 'Value',0, 'Position',[fig.Position(3)/4, 0, fig.Position(3)/2, fig.Position(4)/16], 'Callback',@sliderCallback);
[~, r_max_pos] = max(vecnorm(r_vec_traj, 2, 2));

function sliderCallback(source, ~)

    global r_func f_func g_func f_d_func g_d_func x0 y0 r0_vec v0_vec r_vec_traj r_max_pos

    dtheta = source.Value * pi/180;
    r = r_func(dtheta);
    f = f_func(dtheta);
    g = g_func(dtheta);
    f_d = f_d_func(dtheta);
    g_d = g_d_func(dtheta);
    
    r_vec = f*r0_vec + g*v0_vec;
    v_vec = f_d*r0_vec + g_d*v0_vec;
    
    quiver(0, 0, r_vec(1), r_vec(2), 'g', 'LineWidth',2);
    hold on; grid on;
    quiver(r_vec(1), r_vec(2), v_vec(1), v_vec(2), .1e3, 'b', 'LineWidth',2);
    plot(r_vec_traj(1,:), r_vec_traj(2,:), 'r--');
    plot(x0, y0, 'bo', 'LineWidth', 2);

    plot([0, r_vec_traj(1,r_max_pos)], [0, r_vec_traj(2,r_max_pos)]);
    hold off;
    title(strcat("\theta = ", num2str(atan2(r_vec(2), r_vec(1))*180/pi)));

end % end function
