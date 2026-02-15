function [insar, obs, nObs] = loadInsarData(insar, geo, cmap, useExistingInsar, wrp_plot, unwplot, covplot)

% Function to ingest and subsample InSAR data from pre-prepared *.mat file
%
% Usage: [insar, obs, nObs] = loadInsarData(insar, geo, cmap)
% Input Parameters:
%       insar: structure with insar data and related information
%       geo: structure with local coordinates origin and bounding box
%       cmap: colormaps for plotting
%       useExistingInSAR: exist preloadad InSAR
%       wrp_plot: plot wrapped of data input
%       unw_plot: plot unwrapped of data input
%       covplot: plot covariance matrix of the inputs
%
% Output Parameters:
%       insar: structure with added subsampled data vector and radar look
%       parameters
%       obs: coordinates of observation points after subsampling
%       nObs: number of observation points
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


%% Initialise variables
global outputDir  % Set global variables

nPoints = 0;
LonLat = zeros(0,2);

%% Start loading data
% Create a figure before the for loop

% ---- Default plotting arguments ----
if nargin < 5 || isempty(wrp_plot)
    wrp_plot = false;
end

if nargin < 6 || isempty(unwplot)
    unwplot = false;
end

if nargin < 7 || isempty(covplot)
    covplot = false;
end

% ---- Plot only if requested ----
if covplot
    figure;
    hold on;
end

% Define a color map or choose specific colors for each plot
colors = lines(length(insar)); % This will generate a colormap with as many colors as there are insar inputs

for i = 1:length(insar)
    prepname=insar{i}.dataPath;
    prepname=[prepname(1:end-4),'_prep.mat'];
    if strcmpi(useExistingInsar,'n')
        loadedData = load(insar{i}.dataPath); % load *.mat file

        % Apply bounding box and remove data points outside the AOI
        iOutBox = find(loadedData.Lon < geo.boundingBox(1) | loadedData.Lon > geo.boundingBox(3) | loadedData.Lat > geo.boundingBox(2) | loadedData.Lat < geo.boundingBox(4));
        if sum(iOutBox)>0
            loadedData.Phase(iOutBox) = [];
            loadedData.Lat(iOutBox) = [];
            loadedData.Lon(iOutBox) = [];
            loadedData.Heading(iOutBox) = [];
            loadedData.Inc(iOutBox) = [];
        end

        convertedPhase = (loadedData.Phase / (4*pi) )  * insar{i}.wavelength;    % Convert phase from radians to m
        los = single(-convertedPhase);  % Convert to Line-of-sigth displacement in m
        if sum(isnan(los))>0
            error(['There are NaNs in ',insar{i}.dataPath])
        end
        ll = [single(loadedData.Lon) single(loadedData.Lat)];   % Create Longitude and Latitude 2Xn matrix
        xy = llh2local(ll', geo.referencePoint);    % Transform from geografic to local coordinates

        nPointsThis = size(ll, 1);   % Calculate length of current InSAR data vector
        xy = double([(1:nPointsThis)', xy'*1000]);   % Add ID number column to xy matrix with local coordinates

        % Extract filename to be included in figure names
        [path, name, ext] = fileparts(insar{i}.dataPath);

        %% Plotting data input
        figDir = fullfile(outputDir, 'Figures');
        if ~exist(figDir, 'dir'); mkdir(figDir); end
        % Wrapped
        if wrp_plot
            f = figure('Position',[1 1 700 700], 'Visible','off'); % set 'on' if you want to see it
            plotInsarWrapped(xy, los, insar{i}.wavelength, cmap, name);
            exportgraphics(f, fullfile(figDir, ['Wrapped_' name '.png']), 'Resolution', 300);
            close(f)
        end

        % Unwrapped
        if unwplot
            f = figure('Position',[1 1 700 700], 'Visible','off');
            plotInsarUnwrapped(xy, los, cmap, name);
            exportgraphics(f, fullfile(figDir, ['Unwrapped_' name '.png']), 'Resolution', 300);
            close(f)
        end

        %% Run data vector subsampling using Quadtree and display
        if insar{i}.quadtreeThresh > 0 && length(los) > 20000 %TODO make here more robust, so data size is bigger than 20000 is not downsampled yet..
            [nb, err, nPts, centers, dLos, polys, xLims, yLims] = quadtree(xy, los', insar{i}.quadtreeThresh, 1000, 1); % Run Quadtree on los vector
            c = max(abs([min(los), max(los)])); % Calculate maximu value for symmetric colormap
            caxis([-c c])
            axis equal; axis tight;
            cbar = colorbar; ylabel(cbar, 'Line-of-sight displacement m','FontSize', 14);
            colormap(cmap.redToBlue)
            xlabel('X distance from local origin (m)','FontSize', 14)
            ylabel('Y distance from local origin (m)','FontSize', 14)
            t = title(['Subsampled data. Number of data points used:', num2str(nb)],'FontSize', 18);
            set(t,'Position',get(t,'Position')+[0 1000 0]);
            drawnow
            saveas(gcf, [outputDir,'/Figures/Subsampled_', name, '.png'])



            %% Extract radar look vector information for subsampled points
            disp 'Extracting radar look vector parameters ...'
            [dHeading] = quadtreeExtract(xy, loadedData.Heading, xLims, yLims);    % Extract heading angle values based on quadtree partition
            disp(['Mean heading angle: ', num2str(mean(dHeading)), ' degrees'])
            [dIncidence] = quadtreeExtract(xy, loadedData.Inc, xLims, yLims);  % Extract values based on quadtree partition
            disp(['Mean incidence angle: ', num2str(mean(dIncidence)), ' degrees'])
            disp(['Max and min LoS displacement in m:', num2str(max(dLos)), '  ', num2str(min(dLos))])
            disp 'Extracting observation points height ...'

            pts.xy = [(1:nb)', centers];    % Generate Nx3 [# x y] matrix with local coordinates of data points (centers of Quadtree cells)


            % Include Quadtree results into insar structure
            insar{i}.obs = pts.xy(:,2:3);
            insar{i}.dLos = dLos;
            insar{i}.dHeading = dHeading;
            insar{i}.dIncidence = dIncidence;

        else
            insar{i}.obs = xy(:,2:3);
            insar{i}.dLos = los';
            insar{i}.dHeading = loadedData.Heading';
            insar{i}.dIncidence = loadedData.Inc';
            nb=size(xy,1);
        end

        insarprep=insar{i};
        save(prepname,'-struct','insarprep')
        clear insarprep
    else
        insar{i}=load(prepname);
        nb=size(insar{i}.obs,1);
    end
    
    pts.LonLat = local2llh(insar{i}.obs'/1000, geo.referencePoint)'; % Convert x and y coordinates into Lon Lat coordinates
    
    if covplot    
        %% Create inverse of covariance matrix
        obs = insar{i}.obs;
        
        [X1,X2] = meshgrid(obs(:,1)); % Create square matrices of Xs
        [Y1,Y2] = meshgrid(obs(:,2)); % Create square matrices of Ys
        H = sqrt((X1-X2).^2 + (Y1 - Y2).^2); % Calculate distance between points
        
        % Assign default values if sill, range and nugget are not provided
        if  ~isfield(insar{i},'sillExp');
            disp 'sillExp value not found, assigning default value'
            if contains(insar{i}.dataPath, 'rng')
                insar{i}.sillExp = 0.005;
            elseif contains(insar{i}.dataPath, 'azi')
                insar{i}.sillExp = 0.060;
            elseif contains(insar{i}.dataPath, 'boi')
                insar{i}.sillExp = 0.004;
            end
            
        end
        
        if ~isfield(insar{i},'nugget');
            disp 'nugget value not found, assigning default value'
            if contains(insar{i}.dataPath, 'rng')
                insar{i}.nugget = 0.04;
            elseif contains(insar{i}.dataPath, 'azi')
                insar{i}.nugget = 0.70;
            elseif contains(insar{i}.dataPath, 'boi')
                insar{i}.nugget = 0.04;
            end
    
        end
        
        if ~isfield(insar{i},'range');
            disp 'range value not found, assigning default value'
            if contains(insar{i}.dataPath, 'rng')
                insar{i}.range = 50000;
            elseif contains(insar{i}.dataPath, 'azi')
                insar{i}.range = 50000;
            elseif contains(insar{i}.dataPath, 'boi')
                insar{i}.range = 50000;
            end
            
        end
        
        covarianceMatrix = insar{i}.sillExp * exp(-H/insar{i}.range) + insar{i}.nugget*eye(nb); % Calculate covariance matrix for exponential model with nugget
        % covarianceMatrix_slipberi=insar{i}.sillExp * exp((-3*H)/insar{i}.range)+insar{i}.nugget*eye(nb);
        insar{i}.rcond =rcond(covarianceMatrix);
        insar{i}.invCov = inv(covarianceMatrix); % Calculate inverse of covariance matrix
           
        %%%%%grapth of covariance agains H.
        H_values = H(2:end,1);  % Flatten the H matrix to a vector
        cov_values = covarianceMatrix(2:end,1);  % Flatten the covariance matrix to a vector
        
        colors = lines(9);
        
        % Determine the marker based on the file type
        if contains(insar{i}.dataPath, 'rng')
            markerType = 'x'; % Square for rng.mat
        elseif contains(insar{i}.dataPath, 'azi')
            markerType = 'o'; % Diamond for azi.mat
        elseif contains(insar{i}.dataPath, 'boi')
            markerType = '^'; % Circle for boi.mat
        else
            markerType = 'o'; % Default marker type (circle) if none matches
        end
        
        % Plot the relationship between H and the covariance matrix
        scatter(H_values, cov_values, 'filled', 'Marker', markerType, 'MarkerEdgeColor', colors(i, :)); % Assign a unique color and marker to each dataset
      
        clear X1 X2 Y1 Y2 obs H covarianceMatrix
        insar{i}.ix = nPoints+1:nPoints+nb; % Extract index of data points in obs vector for this interferogram
        nPoints = nPoints + nb; % Number of data points
        LonLat = [LonLat; pts.LonLat]; % Longitude Latitude matrix
    end
end

if covplot
    % Finalize the plot
    xlabel('H (Distance between Points) (m)');
    ylabel('Covariance Matrix Elements');
    %title('Covariance Matrix Elements vs. H for All InSAR Inputs');
    legend(arrayfun(@(i) sprintf('InSAR Input %d', i), 1:length(insar), 'UniformOutput', false));
    grid on;
    hold off; % Release the figure
    xlim([0, 1*10^5]);
    % Save the figure
    saveas(gcf, fullfile(outputDir, 'Figures', 'Covariance_Matrix_vs_H.png')); 
end

obs = llh2local([LonLat'; zeros(1,nPoints)],geo.referencePoint')*1000;
obs = [obs; zeros(1, size(obs,2))]; % Coordinates of observation points
nObs = size(obs,2); % Total number of observation points
