function [theta] = local_sidereal_time(y,m,d,lambda,UT)

J0 = 367*y - double(int64(7/4*(y+int64((m+9)/12)))) + double(int64(275*m/9)) + d + 1721013.5;
T0 = (J0 - 2451545) / 36525;
theta_G0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - 2.583e-8*T0^3;
theta_G0 = theta_G0 - double(int64(theta_G0/360))*360;
theta_G = theta_G0 + 360.98564724*UT/24;
theta_G = theta_G - double(int64(theta_G/360))*360;
theta = theta_G + lambda;
theta = theta - double(int64(theta/360))*360;
end % endfunction
