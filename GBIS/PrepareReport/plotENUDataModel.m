function plotENUDataModel(enu, geo, invpar, invResults, modelinput, saveName, saveflag)

% Function to generate plot with comparison between enu data and model
%
% Usage: plotenuDataModel(enu, geo, invpar, invResults, modelinput, saveName, saveflag)
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
%%
global outputDir  % Set global variables
for i=1:length(enu)
    obsenu = llh2local([enu{i}.ll';zeros(1,length(enu{i}.ll(:,1)))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
    obsenu = [obsenu; zeros(1,size(obsenu,2))]; % Add zeros to third column of observation matrix
    nObsenu = size(obsenu,2); % Determine number of entries in enu observation matrix
    %% Calculate forward model using optimal source parameters

    %obsenu(:,end) = []; % remove coordinates of legend
    %enu.displacements(:,end) = []; % remove displacements for legend
    modenu = forwardGPSModel(obsenu',invpar,invResults,modelinput,i); % Calculate modelled displacements   
    
    % plot E N U observation
    figDir = fullfile(outputDir, 'Figures');
    if ~exist(figDir, 'dir'); mkdir(figDir); end

    % Build xy exactly like InSAR format
    nENU = size(obsenu,2);
    xy = [(1:nENU)' , obsenu(1,:)' , obsenu(2,:)'];

    cmap.redToBlue = crameri('vik');

    compNames = {'E','N','U'};
    enu3 = enu{1}.displacements;   % 3 x N

    % --- one figure, 1x3 layout
    f = figure('Position',[1 1 1600 550], 'Visible','off');
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');
    
    for k = 1:3
        nexttile;
        name_k = [compNames{k}];
        los = enu3(k,:).';   % Nx1
    
        plotGPSscatter(xy, los, cmap, name_k);   % your function should plot into current axes
    end
    
    % save once
    exportgraphics(f, fullfile(figDir, 'obs_E_N_U.png'), 'Resolution', 300);
    close(f);

    %%Plotting Model
    % --- one figure, 1x3 layout
    f = figure('Position',[1 1 1600 550], 'Visible','off');
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');
    
    for k = 1:3
        nexttile;
        name_k = [compNames{k}];
        los = modenu(k,:).';   % Nx1
    
        plotGPSscatter(xy, los, cmap, name_k);   % your function should plot into current axes
    end
    
    % save once
    exportgraphics(f, fullfile(figDir, 'model_E_N_U.png'), 'Resolution', 300);
    close(f);

    %%Plotting Residual
    % --- one figure, 1x3 layout
    f = figure('Position',[1 1 1600 550], 'Visible','off');
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');
    
    for k = 1:3
        nexttile;
        name_k = [compNames{k}];
        residual = enu3(k,:).' - modenu(k,:).';   % observed - model
        los = residual;                           % Nx1
    
        plotGPSscatter(xy, los, cmap, name_k);   % your function should plot into current axes
    end
    
    % save once
    exportgraphics(f, fullfile(figDir, 'residual_E_N_U.png'), 'Resolution', 300);
    close(f);
end
%drawnow

