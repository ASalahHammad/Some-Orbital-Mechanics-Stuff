%% This is algorithm 4.2
function [a,e,i,OMEGA,omega,theta,h] = orbital_elements(r_vec, v_vec)

Re = 6378;
mu = 398600;

r = norm(r_vec);
v = norm(v_vec);

vr = dot(v_vec, r_vec) / r;
h_vec = cross(r_vec, v_vec);
h = norm(h_vec);
i = acosd(h_vec(3)/h);
N_vec = cross([0;0;1], h_vec);
N = norm(N_vec);
if(N_vec(2)>=0)
    OMEGA = acosd(N_vec(1)/N);
else
    OMEGA = 360 - acosd(N_vec(1)/N);
end
e_vec = 1/mu*((v^2-mu/r)*r_vec - r*vr*v_vec);
e = norm(e_vec);
if(e_vec(3)>=0)
    omega = acosd(dot(N_vec,e_vec)/N/e);
else
    omega = 360 - acosd(dot(N_vec,e_vec)/N/e);
end
if(vr>=0)
    theta = acosd(dot(e_vec,r_vec)/e/r);
else
    theta = 360-acosd(dot(e_vec,r_vec)/e/r);
end

rp = h^2/mu/(1+e);
ra = h^2/mu/(1-e);
a = (rp+ra)/2; % this is semimajor axis, not elevation!

end