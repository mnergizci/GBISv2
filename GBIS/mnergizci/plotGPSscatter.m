function plotGpsScatter(xy, vals, cmap, name)
% xy: Nx3 [ID X Y] (m)
% vals: Nx1 (m)
% cmap: struct with field redToBlue (Nx3)

    scatter(xy(:,2), xy(:,3), 10, vals, 'filled'); % 80 = marker size
    colormap(cmap.redToBlue);
    c = max(abs([min(vals), max(vals)]));
    caxis([-c c]);
    axis equal; axis tight;

    ax = gca;
    grid on
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.LineWidth = 1.0;
    ax.GridLineStyle = '--';

    cbar = colorbar; ylabel(cbar,'Displacement (m)','FontSize', 14);
    xlabel('X distance from origin (m)','FontSize', 14)
    ylabel('Y distance from origin (m)','FontSize', 14)
    title(name,'FontSize', 18);
end
