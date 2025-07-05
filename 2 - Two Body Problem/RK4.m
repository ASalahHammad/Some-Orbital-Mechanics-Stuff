function [t, y] = RK4(f_dot, t, y0)

y = nan(length(y0), length(t));
y(:, 1) = y0;
for n = 1:length(t)-1
    h = t(n+1) - t(n);
    K1 = h*f_dot(t(n), y(:, n));
    K2 = h*f_dot(t(n)+h/2, y(:, n)+0.5*K1);
    K3 = h*f_dot(t(n)+h/2, y(:, n)+0.5*K2);
    K4 = h*f_dot(t(n)+h, y(:, n)+K3);
    y(:, n+1) = y(:, n) + 1/6*(K1 + 2*K2 + 2*K3 + K4);
end

y = y.';
end % endfunction
