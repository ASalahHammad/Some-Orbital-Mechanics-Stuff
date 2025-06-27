%% Algorithm 5.1
clear
close all
clc

%% All units in km
mu = 398600;
Re = 6378;

r1_vec = [-294.32; 4265.1; 5986.7];
r2_vec = [-1365.5; 3637.6; 6346.8];
r3_vec = [-2940.3; 2473.7; 6555.8];

r1 = norm(r1_vec);
r2 = norm(r2_vec);
r3 = norm(r3_vec);
ur1 = r1_vec/r1;

C12 = cross(r1_vec, r2_vec);
C23 = cross(r2_vec, r3_vec);
C31 = cross(r3_vec, r1_vec);
C23_hat = C23/norm(C23);

%% check
if (abs(dot(ur1,C23_hat)) <=1e-3)
    fprintf("Three vectors are coplanar\n");
    N = r1*C23 + r2*C31 + r3*C12;
    D = C12 + C23 + C31;
    S = r1_vec*(r2-r3) + r2_vec*(r3-r1) + r3_vec*(r1-r2);
    v2_vec = sqrt(mu/norm(N)/norm(D)) * (cross(D,r2_vec)/r2 + S);
    v2 = norm(v2_vec);
    %% Use Algorithm 4.2
    [a,e,i,OMEGA,omega,theta] = orbital_elements(r2_vec,v2_vec)
    h = norm(cross(r2_vec, v2_vec));
    rp = h^2/mu/(1+e);
    T = 2*pi/sqrt(mu)*a^(3/2) / 3600; % hr
    % perigee altitude
    perigee_alt = rp - Re;
    % time since perigee
    E1 = 2*atan(sqrt((1-e)/(1+e))*tand(theta/2)); % rad
    Me1 = E1 - e*sin(E1);
    t1 = h^3/mu^2/(1-e^2)^(3/2)*Me1;
end
