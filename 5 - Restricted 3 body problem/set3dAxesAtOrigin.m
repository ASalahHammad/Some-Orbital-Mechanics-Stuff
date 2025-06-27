function [qq] = set3dAxesAtOrigin(ax)

xrange = xlim;
yrange = ylim;
zrange = zlim;

if(nargin>0)
    qq = quiver3(ax, [0,0,0], [0,0,0], [0,0,0], [xrange(2)-xrange(1),0,0], [0,yrange(2)-yrange(1),0], [0,0,zrange(2)-zrange(1)], 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % X-axis
else
    qq = quiver3([0,0,0], [0,0,0], [0,0,0], [xrange(2)-xrange(1),0,0], [0,yrange(2)-yrange(1),0], [0,0,zrange(2)-zrange(1)], 'k', 'LineWidth', 1.5, 'AutoScale', 'off'); % X-axis
end

end
