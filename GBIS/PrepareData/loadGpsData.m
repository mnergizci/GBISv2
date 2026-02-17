function [gps, obsGps, nObsGps] = loadGpsData(gps, geo, gpsplot, enuplot)

% Function to ingest GPS data from text files
%
% Usage: [gps, obsGps, nObsGps] = loadGpsData(gps, geo)
% Input Parameters:
%       gps: gps structure containing path to data files
%       geo: structure with local coordinates origin and bounding box
%
% Output Parameters:
%       gps: structure with GPS data and further information (e.g., inverse
%       of covariance matrix)
%       obsGps: coordinates of observation points
%       nObsGps: number of observation points
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

% ---- Default plotting arguments ----
if nargin < 3 || isempty(gpsplot)
    gpsplot = false;
end

if nargin < 4 || isempty(enuplot)
    enuplot = false;
end

global outputDir  % Set global variables

gpsTxt = load(gps.dataPath); % Read text file with GPS data

gps.ll = gpsTxt(:, 1:2); % Read GPS site Longitude and Latitude coordinates 
nGps = size(gpsTxt, 1); % Retrieve number of GPS sites
gps.displacements = gpsTxt(:, [5,3,7])'/1000; % Import GPS displacements in mm and tranforms to m
gps.sigmas = gpsTxt(:, [6,4,8])'/1000; % Import GPS displacement st. dev. in mm and tranforms to m
gps.variance = gps.sigmas.^2; % Calculate variance of GPS displacements
gps.invCov = diag(1./reshape(gps.variance(1:3,:),nGps*3,1)); % Generate inverse of covariance matrix for GPS

obsGps = llh2local([gps.ll'; zeros(1,nGps)], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
obsGps = [obsGps; zeros(1,size(obsGps,2))]; % Add zeros to third column of observation matrix
nObsGps = size(obsGps,2); % Determine number of entries in GPS observation matrix

if gpsplot
    % Display Gps vectors
    obsGps(:,end+1) = [max(obsGps(1,:))+5000; min(obsGps(2,:))-5000; 0]; % add coordinates of legend
    scalebar = abs(round(max(gps.displacements(:))/3,3));
    gps.displacements(:,end+1) = [-scalebar 0 0]; % add displacements for legend
    
    figure
    quiver(obsGps(1,:), obsGps(2,:), gps.displacements(1,:), gps.displacements(2,:), 1, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 0.03, 'Marker', 's')
    axis equal; 
    ax = gca;
    grid on
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.LineWidth = 1.5;
    ax.GridLineStyle = '--';
    xlabel('X distance from local origin (m)')
    ylabel('Y distance from local origin (m)')
    title('GPS horizontal displacements')
    xlim([min(obsGps(1,:))-10000 max(obsGps(1,:))+10000]);
    ylim([min(obsGps(2,:))-10000 max(obsGps(2,:))+10000]);
    text(obsGps(1,end),obsGps(2,end)-2000,[num2str(scalebar*1000),' mm'])
    drawnow
    saveas(gcf,[outputDir,'/Figures/GPS_displacements.png'])
    obsGps(:,end) = []; % remove coordinates of legend
    gps.displacements(:,end) = []; % remove displacements for legend
end

if enuplot

    figDir = fullfile(outputDir, 'Figures');
    if ~exist(figDir, 'dir'); mkdir(figDir); end

    % Build xy exactly like InSAR format
    nGps = size(obsGps,2);
    xy = [(1:nGps)' , obsGps(1,:)' , obsGps(2,:)'];

    cmap.redToBlue = crameri('vik');

    compNames = {'E','N','U'};
    enu = gps.displacements;   % 3 x N

    for k = 1:3
        name = ['GPS_' compNames{k}];

        f = figure('Position',[1 1 700 700], 'Visible','off');

        los = enu(k,:).';   % MUST be Nx1
        plotGPSscatter(xy, los, cmap, name);

        exportgraphics(f, fullfile(figDir, [name '.png']), 'Resolution', 300);
        close(f);
    end
end








