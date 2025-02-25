function [] = fitVariogram(inputFile, wavelength, applyfiltering, pos_path)
% Function to fit an exponential function to the data isotropic (semi-)variogram
%
% Usage:  fitVariogram(inputFile, wavelength)
%   inputFile:  path and name of *.mat file containing data in the format used by
%               GBIS
%
%   wavelength: wavelength of InSAR data in meters (e.g., 0.056 m for
%               Sentine-1/Envisat/ERS)
%
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
close all
clc

if nargin <4
    pos_path = [];
end
if nargin < 3
    applyfiltering = 'No';
end 

% Load dataset
disp('Ingesting data to estimate (semi-)variogram ...')
insarData = load(inputFile);

% Create colormap
cmapSeismo = colormap_cpt('GMT_seis.cpt', 256);    % GMT 'Seismo' colormap for wrapped data

% % Create colormap
% cmapRoma = colormap_cpt('GMT_roma.cpt', 100);    % GMT 'Roma' colormap for wrapped data
% 

% Find a local reference point
refPoint = [min(insarData.Lon), min(insarData.Lat)]; % Determine local reference point

% Convert phase to LOS displacement
convertedPhase = (insarData.Phase / (4*pi)) * wavelength;   % Convert phase from radians to m
los = single(-convertedPhase);                              % Convert to Line-of-sigth displacement in m

% Determine subsampling factor for faster plotting
if length(los) > 400000 && length(los) < 1000000
    sampling = 2;
elseif length(los) > 1000000
    sampling = 5;
else
    sampling = 1;
end


% % Plot wrapped dataset
% figure
% scatter(insarData.Lon(1:sampling:end), insarData.Lat(1:sampling:end), [], mod(los(1:sampling:end), wavelength/2),'.');
% colormap(cmapSeismo)
% caxis([0 wavelength/2])
% axis xy
% axis equal
% title('Wrapped Interferogram')
% xlabel('Longitude (degrees)')
% ylabel('Latitude (degrees)')
% colorbar

% This is much better for me?
figure;
scatter(insarData.Lon(1:sampling:end), insarData.Lat(1:sampling:end), [], los(1:sampling:end), '.');
colormap(cmapSeismo);
caxis([-1 1]);
axis xy;
axis equal;
title('Unwrapped Interferogram');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
colorbar;


% Debugging: check if pos_path is empty or not
if ~isempty(pos_path)
    disp('Using provided pos_path.');
else
    disp('pos_path is empty, entering else section.');
end

%% Give choice if semi-variogram should be calculated in rectangular area of interest or over the entire image with a mask
%% Use the provided polygon coordinates if given, otherwise prompt for user input
if ~isempty(pos_path)
    % Use the provided polygon coordinates
    load(pos_path);
    if exist('pos', 'var') && ~isempty(pos)
        disp('Using provided polygon coordinates.');
        in = inpolygon(insarData.Lon, insarData.Lat, pos(:,1), pos(:,2));
        ixSubset = find(in == 0); % Points outside the polygon

    elseif exist('pos_rectangular', 'var') && ~isempty(pos_rectangular)
        disp('Using provided rectangular coordinates.');
        % Extract the min and max values from the saved rectangular coordinates
        maxLon = pos_rectangular(1);
        maxLat = pos_rectangular(2);
        minLon = pos_rectangular(3);
        minLat = pos_rectangular(4);
        ixSubset = find(insarData.Lat > minLat & insarData.Lat < maxLat & insarData.Lon < minLon & insarData.Lon > maxLon);
    else
        error('The provided file does not contain valid polygon or rectangular coordinates.');
    end

else
    % Prompt user to select region of interest
    choice = questdlg('Would you like to select a rectangular area or mask out a region?', 'Option:', 'Rectangle', 'Mask','Rectangle');

    switch choice
        case 'Rectangle'
            disp('Select rectangular area using mouse.');
            bounds = getrect;
            maxLon = bounds(1);
            minLon = bounds(1)+bounds(3);
            minLat = bounds(2);
            maxLat = bounds(2)+bounds(4);
            pos_rectangular = [maxLon; maxLat; minLon; minLat];
            save('pos_rectangular.mat', 'pos_rectangular');
            ixSubset = find(insarData.Lat > minLat & insarData.Lat < maxLat & insarData.Lon < minLon & insarData.Lon > maxLon);
        case 'Mask'
            disp('Draw closed polygon using mouse.');
            polyMask = impoly;
            pos = getPosition(polyMask);
            save('pos.mat', 'pos');
            in = inpolygon(insarData.Lon, insarData.Lat, pos(:,1), pos(:,2));
            ixSubset = find(in == 0);
    end
end


% Extract subset from rectangular area or after masking
subset = los(ixSubset);
llon = insarData.Lon(ixSubset);
llat = insarData.Lat(ixSubset);

% Display subregion from selection
figure('Position', [1, 1, 1200, 1000]);
subplot(2,3,1)
% scatter(llon(:),llat(:),[],mod(subset(:),wavelength/2),'.')
% colormap(cmapSeismo)
% caxis([0 wavelength/2])
scatter(llon(:), llat(:), [], subset(:), '.');
colormap(cmapSeismo);
caxis([-1 1]);


axis xy
axis equal
axis tight
title('NON-DETRENDED')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

%% Remove linear trend from subregion
sll = [llon'; llat']';
xy = llh2local(sll',refPoint);
xy = xy*1000;

A = [xy' ones([length(xy) 1])];

coeff = lscov(A,subset);
deramped = subset - A*coeff;

%% Display trend and subregion after removal of trend
subplot(2,3,2)
% scatter(llon(:),llat(:),[],mod(A(:,:)*coeff,wavelength/2),'.')
% colormap(cmapSeismo)
% caxis([0 wavelength/2])
scatter(llon(:), llat(:), [], A(:,:) * coeff, '.');
colormap(cmapSeismo);
caxis([-0.5 0.5]);
axis xy
axis equal
axis tight
title('ESTIMATED TREND')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

if strcmp(applyfiltering, 'Yes')
    % Calculate mean and standard deviation
    mu = mean(deramped);
    sigma = std(deramped);

    % Identify outliers using the 3-sigma rule
    outliers = deramped > (mu +sigma) | deramped < (mu -sigma);

    % Remove outliers based on the 3-sigma rule
    deramped = deramped(~outliers);
    llon_new = llon(~outliers);
    llat_new = llat(~outliers);
    sll_new = [llon_new'; llat_new']';
    xy_new = llh2local(sll_new',refPoint);
    xy_new = xy_new*1000;
    
end

subplot(2,3,3)
% scatter(llon(:),llat(:),[],mod(deramped(:),wavelength/2),'.')
% colormap(cmapSeismo)
% caxis([0 wavelength/2])

scatter(llon_new(:),llat_new(:),[],deramped(:),'.')
colormap(cmapSeismo);
caxis([-1 1]);
axis xy
axis equal
axis tight
title('DETRENDED')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

%% Calculate and display variogram before plane removal
subplot(2,3,4)
variog = variogram(xy',double(subset),'plotit',true,'subsample',10000,'nrbins',40);
title('NON-DETRENDED')

% Calculate and display variogram after detrending
variogDtrnd = variogram(xy_new',double(deramped),'plotit',false,'subsample',10000,'nrbins',40);

%% Fit exponential function to experimental variogram and display
subplot(2,3,5)
[a,c,n] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 3);
title('DETRENDED')
disp('done')

h =subplot(2,3,6);
set(h,'visible','off')
text(0.1,1.0,'Fitted exponential semi-variogram parameters:','FontSize',14)
text(0.1,0.8,['Sill (minus nugget):  ', num2str(c)],'FontSize',14)
text(0.1,0.6,['Range:  ', num2str(a)],'FontSize',14)
text(0.1,0.4,['Nugget:  ', num2str(n)],'FontSize',14)

% Print variogram exponential fit parameters to screen
disp(['Sill:  ',num2str(c)])
disp(['Range:  ',num2str(a)])
disp(['Nugget:  ',num2str(n)])

% Extract the base name of the .mat file and create the .cov filename
[~, baseName, ~] = fileparts(inputFile);
covFileName = [baseName, '.cov'];

% Save the main figure
saveas(gcf, [baseName, '_semivari_all_together.png']);  % Saves the figure as a PNG file

% Write the parameters to the .cov file
fileID = fopen(covFileName, 'w');
fprintf(fileID, 'Sill: %f\n', c);
fprintf(fileID, 'Range: %f\n', a);
fprintf(fileID, 'Nugget: %f\n', n);
fclose(fileID);

disp(['Parameters written to ', covFileName]);



%%%%extra_figure:
extraFigure = figure;
[a,c,n] = variogramfitalone(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 3);
% Save the extra figure
saveas(extraFigure, [baseName, '_semivariogram_noredline.png']);  % Saves the extra figure as a PNG file