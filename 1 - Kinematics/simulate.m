function [] = simulate(t, r, v, a)

r_c = nan(size(r, 1), length(t));

figure; hold on; grid on;
xlabel("x"); ylabel("y"); zlabel("z");
view([45, 45]);
% min_xyz = -100;
% max_xyz = 100;
% xlim([min_xyz, max_xyz]);
% ylim([min_xyz, max_xyz]);
% zlim([min_xyz, max_xyz]);
axis equal

pp1 = plot(nan, nan);
pp2 = plot(nan, nan);
pp3 = quiver3(nan, nan, nan, nan, nan, nan);
pp4 = quiver3(nan, nan, nan, nan, nan, nan);
pp5 = quiver3(nan, nan, nan, nan, nan, nan);
pp6 = quiver3(nan, nan, nan, nan, nan, nan);
pp7 = quiver3(nan, nan, nan, nan, nan, nan);
pp8 = quiver3(nan, nan, nan, nan, nan, nan);
pp9 = quiver3(nan, nan, nan, nan, nan, nan);
pp10 = quiver3(nan, nan, nan, nan, nan, nan);

for n = 1:length(t)
    if(exist("a", "var"))
        u_t = v(:, n)/norm(v(:, n));
        u_b = cross(v(:, n), a(:, n))/norm(cross(v(:, n), a(:, n)));
        u_n = cross(u_b, u_t);
        a_n = dot(a(:, n), u_n);
        rho = dot(v(:, n), v(:, n))/a_n;
        r_c(:, n) = r(:, n) + rho*u_n;
    end

    % Plot
    scale_factor = 8;

    delete(pp1); delete(pp2); delete(pp3); delete(pp4); delete(pp5); delete(pp6); delete(pp7); delete(pp8); delete(pp9); delete(pp10);
    pp1 = quiver3(0, 0, 0, r(1, n), r(2, n), r(3, n), 'g', 'LineWidth', 0.5, 'AutoScale', 'off');
    pp2 = plot3(r(1, 1:n), r(2, 1:n), r(3, 1:n), 'b', 'LineWidth', 2);
    pp3 = quiver3(r(1, n), r(2, n), r(3, n), u_n(1)*scale_factor, u_n(2)*scale_factor, u_n(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
    pp4 = quiver3(r(1, n), r(2, n), r(3, n), u_t(1)*scale_factor, u_t(2)*scale_factor, u_t(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
    pp5 = quiver3(r(1, n), r(2, n), r(3, n), u_b(1)*scale_factor, u_b(2)*scale_factor, u_b(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
    pp6 = plot3(r(1, n), r(2, n), r(3, n), 'bo', 'LineWidth', 2);
    if(exist("a", "var"))
        pp7 = plot3(r_c(1, 1:n), r_c(2, 1:n), r_c(3, 1:n), 'ro', 'LineWidth', 2);
        pp8 = quiver3(r(1, n), r(2, n), r(3, n), r_c(1, n)-r(1, n), r_c(2, n)-r(2, n), r_c(3, n)-r(3, n), 'k', 'LineWidth', 0.5, 'AutoScale', 'off');
        % pp9 = quiver3(0, 0, 0, r_c(1, n), r_c(2, n), r_c(3, n), 'c', 'LineWidth', 0.5);
        if(n~=1)
            pp10 = fill3([r(1, n), r_c(1, n), r_c(1, n-1), r(1, n-1)], [r(2, n), r_c(2, n), r_c(2, n-1), r(2, n-1)], [r(3, n), r_c(3, n), r_c(3, n-1), r(3, n-1)], 'm');
        end
    end
    title(strcat("Time: ", num2str(t(n)), " sec"), 'Interpreter', 'latex');
    set3dAxesAtOrigin();
    pause(0.02);
end
legend("Position Vector", "Trajectory", "unit normal", "unit tangent", "unit binormal", "R_c");

end % endfunction
