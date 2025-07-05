function [] = simulate(t, dt, r, v, a, styles, legends, ax, obj)

N_bodies = int32(size(r, 1) / 3);

% dt is the jump in the tme stamps, it's not the actual delta_time
if(~exist("styles", "var"))
    styles = ["r", "g", "b"];
end
if(~exist("legends", "var"))
    legends = compose("body %i", 1:N_bodies);
end
if(~exist("ax", "var"))
    fig = figure;
    ax = gca;
end
if(exist("obj", "var"))
    for k = 1:length(obj)
        obj(k).objImage = imread(obj(k).filename);
        [X_im, Y_im, Z_im] = sphere;
        obj(k).X_im = obj(k).radius * X_im;
        obj(k).Y_im = obj(k).radius * Y_im;
        obj(k).Z_im = obj(k).radius * Z_im;
    end
end

body(N_bodies) = struct("colour", "b");

for i = 1:N_bodies
    body(i).style = styles(i); % colours and styles
end

r_c = nan(size(r, 1), length(t));

axes(ax);
hold on; grid on;
xlabel("x"); ylabel("y"); zlabel("z");
view([150, 20]);


axis equal
% pbaspect([1,1,1]);
pp1 = plot3(ax, nan, nan, nan, nan, nan, nan);
pp2 = plot3(ax, nan, nan, nan, nan, nan, nan);
pp3 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp4 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp5 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp6 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp7 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp8 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp9 = quiver3(ax, nan, nan, nan, nan, nan, nan);
pp10 = quiver3(ax, nan, nan, nan, nan, nan, nan);
qq = quiver3(ax, nan, nan, nan, nan, nan, nan);

for n = 1:dt:length(t)
    delete(qq); delete(pp1); delete(pp2); delete(pp3); delete(pp4); delete(pp5); delete(pp6); delete(pp7); delete(pp8); delete(pp9); delete(pp10);
    for k = 1:N_bodies
        ind_x = (k-1)*3+1;
        ind_y = (k-1)*3+2;
        ind_z = (k-1)*3+3;
        if(exist("a", "var"))
            u_t = v(ind_x:ind_z, n)/norm(v(ind_x:ind_z, n));
            u_b = cross(v(ind_x:ind_z, n), a(ind_x:ind_z, n))/norm(cross(v(ind_x:ind_z, n), a(ind_x:ind_z, n)));
            u_n = cross(u_b, u_t);
            a_n = dot(a(ind_x:ind_z, n), u_n);
            rho = dot(v(ind_x:ind_z, n), v(ind_x:ind_z, n))/a_n;
            r_c(ind_x:ind_z, n) = r(ind_x:ind_z, n) + rho*u_n;
        end

        % Plot
        scale_factor = 8;
    
        pp1(k) = quiver3(ax, 0, 0, 0, r(ind_x, n), r(ind_y, n), r(ind_z, n), body(k).style, 'LineWidth', 0.5, 'AutoScale', 'off');
        pp2(k) = plot3(ax, r(ind_x, 1:n), r(ind_y, 1:n), r(ind_z, 1:n), body(k).style, 'LineWidth', 2);
        if(exist("obj", "var") && k<=length(obj))
            pp3(k) = surf(ax, obj(k).X_im+r(ind_x,n), obj(k).Y_im+r(ind_y,n), obj(k).Z_im+r(ind_z,n), 'CData', flipud(obj(k).objImage), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        else
            pp3(k) = plot3(ax, r(ind_x, n), r(ind_y, n), r(ind_z, n), strcat(body(k).style,'o'), 'LineWidth', 2);
        end
        if(exist("a", "var"))
            % pp4(k) = quiver3(ax, r(ind_x, n), r(ind_y, n), r(ind_z, n), u_t(1)*scale_factor, u_t(2)*scale_factor, u_t(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
            % pp5(k) = quiver3(ax, r(ind_x, n), r(ind_y, n), r(ind_z, n), u_b(1)*scale_factor, u_b(2)*scale_factor, u_b(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
            % pp6(k) = quiver3(ax, r(ind_x, n), r(ind_y, n), r(ind_z, n), u_n(1)*scale_factor, u_n(2)*scale_factor, u_n(3)*scale_factor, 'r', 'LineWidth', 1, 'AutoScale', 'off');
            % pp7(k) = plot3(ax, r_c(ind_x, 1:n), r_c(ind_y, 1:n), r_c(ind_z, 1:n), 'ro', 'LineWidth', 0.2); % centre of rotation
            % pp8(k) = quiver3(ax, r(ind_x, n), r(ind_y, n), r(3, n), r_c(ind_x, n)-r(ind_x, n), r_c(ind_y, n)-r(ind_y, n), r_c(ind_z, n)-r(ind_z, n), 'k', 'LineWidth', 0.5, 'AutoScale', 'off'); % arrow connecting r to r_c
            % pp9(k) = quiver3(ax, 0, 0, 0, r_c(ind_x, n), r_c(ind_y, n), r_c(ind_z, n), 'c', 'LineWidth', 0.5, 'AutoScale', 'off'); % position of centre of rotation
            if(n~=1)
                % pp10(k) = fill3([r(ind_x, n), r_c(ind_x, n), r_c(ind_x, n-1), r(ind_x, n-1)], [r(ind_y, n), r_c(ind_y, n), r_c(ind_y, n-1), r(ind_y, n-1)], [r(ind_z, n), r_c(ind_z, n), r_c(ind_z, n-1), r(ind_z, n-1)], 'm'); % plane of rotation
            end
        end

    end
    title(ax, strcat("Time: ", num2str(t(n)), " sec"), 'Interpreter', 'latex');
    qq = set3dAxesAtOrigin(ax);
    pause(0.01);
end
legend(ax, pp3, legends);

end % endfunction
