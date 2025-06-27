clear
close all
clc

function sliderCallback(source, ~)
    e = source.Value;
    m1 = 5.974e24; % earth
    m2 = 641.9e21; % moon
    G = 6.6743e-11;
    mu = G*(m1+m2);

    a = 384.4e6; % semi major axis of moon
    r_a = a*(1+e);
    h = sqrt(mu*r_a);
    % theta = linspace(0, 2*pi, 80);
    theta = 0.5*(1 - cos(linspace(0, pi, 50))) * 2*pi;
    r = h^2 / mu ./ (1+e*cos(theta));
    plot([0, 0], 'bo', 'LineWidth', 3);
    hold on;
    plot(r.*cos(theta), r.*sin(theta), 'r--', 'LineWidth', 2);
    % hold off;
    title(strcat("e = ", num2str(e)));
    pbaspect([1 1 1]);
    axis equal;
    axis off;
    xlim([-1.5e9, 1.5e9]);
    ylim([-9e8, 9e8]);
end

temp.Value = 1;

fig = figure('Position', [1, 1, 1366, 650]);
sliderCallback(temp);
clear temp
slider = uicontrol("Style","slider", "Min",0, "Max", 2, 'Value',1, 'Position',[fig.Position(3)/4, 0, fig.Position(3)/2, fig.Position(4)/16], 'Callback',@sliderCallback);

%% Some other Plots
% Z = 0:100:1e6;
% v = sqrt(mu./(RE + Z));
% T = 2*pi/sqrt(mu)*(RE + Z).^1.5 / 60;
% figure;
% subplot(1,2,1); plot(Z, v); xlabel("Z [m]"); ylabel("v [m/s]"); grid on;
% subplot(1,2,2); plot(Z, T); xlabel("Z [m]"); ylabel("T [min.]"); grid on;
