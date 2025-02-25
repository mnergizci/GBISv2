function plotInsarWrapped(xy, los, wavelength, cmap, name)

% Function to plot wrapped InSAR data
%
% Usage: plotInsarWrapped(xy, los, cmap, name)
% Input Parameters:
%       xy: local coordinates of data points
%       los: line-of-sight displacement measured at data points
%       wavelength: data wavelength for wrapping
%       cmap: colormaps for plotting
%       name: name of dataset for figure title
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Email: gbis.software@gmail.com
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
% Last update: 8 August, 2018

%% Patch scattered data for faster plotting
%     edge = round(min(abs(diff(xy(:,3)))))+2; % Size of patch set to minumum distance between points
%     
%     if edge < 400
%         edge = 400;
%     end
%     
%     xy=xy/1000; % convert to km
%     edge=edge/1000;
% 
%     xs = [xy(:,2)'; xy(:,2)'+edge; xy(:,2)'+edge; xy(:,2)']; % Coordinates of four vertex of patch
%     ys = [xy(:,3)'; xy(:,3)'; xy(:,3)'+edge; xy(:,3)'+edge];
%      
%    
%     % Display wrapped interferogram at 5.6 cm wavelength
%     colormap(cmap.seismo); 
%     h1 = patch(xs, ys, 'r');
%     set(h1, 'facevertexcdata', mod(los,wavelength/2), 'facecolor', 'flat', 'edgecolor', 'none')
%     axis equal; axis tight;
%     ax = gca;
%     grid on
%     ax.Layer = 'top';
%     ax.Box = 'on';
%     ax.LineWidth = 1.0;
%     ax.GridLineStyle = '--';
%     %cbar = colorbar; ylabel(cbar,'Line-of-sight displacement m','FontSize', 14); 
%     xlabel('X distance from local origin (m)','FontSize', 14)
%     ylabel('Y distance from local origin (m)','FontSize', 14)
%     t = title(['Wrapped InSAR Data: ', name],'FontSize', 18);
%     set(t,'Position',get(t,'Position')+[0 1000 0]);
%     drawnow
        % Convert coordinates to km
        xy = xy / 1000; % convert to km

        % Compute mean and standard deviation of line-of-sight (los) displacement
        mean_los = mean(los);
        std_los = std(los);

        % Set colorbar limits to ±4 from the minimum and maximum values
%         cmin = round(min(los) - 4);
%         cmax = round(max(los) + 4);
        cmin = round(-4);
        cmax = round(4);

        % Display unwrapped interferogram
        colormap(cmap.redToBlue);

        % Create scatter plot
        scatter(xy(:,2), xy(:,3), 20, los, 'filled');
        axis equal; axis tight;
        ax = gca;
        grid on;
        ax.Layer = 'top';
        ax.Box = 'on';
        ax.LineWidth = 1.0;
        ax.GridLineStyle = '--';

        % Set color limits and colorbar
        caxis([cmin cmax]);
        cbar = colorbar;
        ylabel(cbar, '', 'FontSize', 14);

        % Set labels and title
%         xlabel('X distance from local origin (km)', 'FontSize', 14);
%         ylabel('Y distance from local origin (km)', 'FontSize', 14);
%         t = title(['Unwrapped InSAR Data: ', name], 'FontSize', 18);
%         set(t, 'Position', get(t, 'Position') + [0 100 0]);

        drawnow;

    

