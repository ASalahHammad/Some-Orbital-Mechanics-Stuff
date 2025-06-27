function [Fdot] = fdot(t, states)

Fdot = nan(size(states));
Fdot(1:3) = states(4:6);
Fdot(4:6) = [16; 30*t; 18/5*t^2+4];

end % endfunction
