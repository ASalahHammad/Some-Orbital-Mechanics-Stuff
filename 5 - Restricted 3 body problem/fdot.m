function [states_dot] = fdot(t, states_0, mass_vec)

G = 6.67e-11; % 

N = length(mass_vec); % no. of bodies, it's different from the number of states

states_0 = reshape(states_0, [3, 2*N]);
R_vec = states_0(:, 1:N);
% Rd_vec = states_0(:, N+1:end);

states_dot = nan(size(states_0));
states_dot(:, 1:N) = states_0(:, N+1:2*N);

for i = 1:N
    array = 1:N; % indices of all bodies except current one
    array = array(array~=i);
    R_i = R_vec(:,i);
    states_dot(:,N+i) = sum( G*mass_vec(array).*(R_vec(:,array) - R_i) ./ vecnorm(R_vec(:,array) - R_i).^3, 2);
end

states_dot = reshape(states_dot, [2*3*N, 1]);

end % endfunction
