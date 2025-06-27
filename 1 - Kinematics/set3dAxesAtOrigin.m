function set3dAxesAtOrigin(ax)

persistent qq1 qq2 qq3
if(exist("qq1", "var"))
    delete(qq1); delete(qq2); delete(qq3);
end
    
xrange = xlim;
yrange = ylim;
zrange = zlim;

if(nargin>0)
    qq1 = quiver3(ax, 0, 0, 0, xrange(2)-xrange(1), 0, 0, 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % X-axis
    qq2 = quiver3(ax, 0, 0, 0, 0, yrange(2)-yrange(1), 0, 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % Y-axis
    qq3 = quiver3(ax, 0, 0, 0, 0, 0, zrange(2)-zrange(1), 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % Z-axis
else
    qq1 = quiver3(0, 0, 0, xrange(2)-xrange(1), 0, 0, 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % X-axis
    qq2 = quiver3(0, 0, 0, 0, yrange(2)-yrange(1), 0, 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % Y-axis
    qq3 = quiver3(0, 0, 0, 0, 0, zrange(2)-zrange(1), 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % Z-axis        
end

end
