function [states_dot] = fdot(t, states_0, m1, m2, G)

if(~exist("G", "var"))
    G = 6.67e-11;
end
X1 = states_0(1);
Y1 = states_0(2);
Z1 = states_0(3);
X2 = states_0(4);
Y2 = states_0(5);
Z2 = states_0(6);

r = sqrt((X2-X1)^2 + (Y2-Y1)^2 + (Z2-Z1)^2);

states_dot = nan(size(states_0));
states_dot(1:6) = states_0(7:12);
states_dot(7) = G*m2*(X2-X1)/r^3;
states_dot(8) = G*m2*(Y2-Y1)/r^3;
states_dot(9) = G*m2*(Z2-Z1)/r^3;
states_dot(10) = -G*m1*(X2-X1)/r^3;
states_dot(11) = -G*m1*(Y2-Y1)/r^3;
states_dot(12) = -G*m1*(Z2-Z1)/r^3;

end % endfunction
