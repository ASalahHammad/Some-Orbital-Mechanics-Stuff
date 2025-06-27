%% Algorithm 4.5
% x_bar means perifocal frame
% X means ECEF frame

function [r_X, v_X] = perifocal2ECEF(h,e,i,OMEGA,omega,theta)

mu = 398600;

r_x_bar = h^2/mu/(1+e*cosd(theta)) * [cosd(theta); sind(theta); 0];
v_x_bar = mu/h * [-sind(theta); e+cosd(theta); 0];

Q_xdash_X = [-sind(OMEGA)*cosd(i)*sind(omega)+cosd(OMEGA)*cosd(omega), -sind(OMEGA)*cosd(i)*cosd(omega)-cosd(OMEGA)*sind(omega), sind(OMEGA)*sind(i);
             cosd(OMEGA)*cosd(i)*sind(omega)+sind(OMEGA)*cosd(omega), cosd(OMEGA)*cosd(i)*cosd(omega)-sind(OMEGA)*sind(omega), -cosd(OMEGA)*sind(i);
             sind(i)*sind(omega), sind(i)*cosd(omega), cosd(i)];
r_X = Q_xdash_X * r_x_bar;
v_X = Q_xdash_X * v_x_bar;


