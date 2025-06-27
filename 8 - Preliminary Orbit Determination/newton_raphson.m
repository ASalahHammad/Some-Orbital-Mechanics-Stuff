function [x2] = newton_raphson(fHandle, fdHandle, x1, tolerance)

if(~exist("tolerance","var"))
    tolerance = 1e-6;
end

epsilon = 1;
while(epsilon >= tolerance)
    f = fHandle(x1);
    fd = fdHandle(x1);
    x2 = x1 - f/fd;
    epsilon = x2-x1;
    x1 = x2;
end

end % endfunction
