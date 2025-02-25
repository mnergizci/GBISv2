function plotGPSDataModel(gps, geo, invpar, invResults, modelinput, saveName, saveflag)

% Function to generate plot with comparison between GPS data and model
%
% Usage: plotGPSDataModel(gps, geo, invpar, invResults, modelinput, saveName, saveflag)
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
for i=1:length(gps)
    obsGPS = llh2local([gps{i}.ll';zeros(1,length(gps{i}.ll(:,1)))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
    obsGPS = [obsGPS; zeros(1,size(obsGPS,2))]; % Add zeros to third column of observation matrix
    nObsGPS = size(obsGPS,2); % Determine number of entries in GPS observation matrix
    %% Calculate forward model using optimal source parameters

    %obsGPS(:,end) = []; % remove coordinates of legend
    %gps.displacements(:,end) = []; % remove displacements for legend
    modGPS = forwardGPSModel(obsGPS',invpar,invResults,modelinput,i); % Calculate modelled displacements   
    
    %% Generate plot
    % Display GPS vectors
    obsGPS(:,end+1) = [min(obsGPS(1,:)); min(obsGPS(2,:))-5000; 0]; % add coordinates of legend
    scalebar = abs(ceil(max(gps{i}.displacements(:))*1000/15))*5/1000; % Determine length of scalebar
    gps{i}.displacements(:,end+1) = [scalebar 0 0]; % add "displacements" of legend

    % Set scale for vectors
    maxobs=max([sqrt(gps{i}.displacements(1,:).^2+gps{i}.displacements(2,:).^2),gps{i}.displacements(3,:)]);
    maxmodel=max([sqrt(modGPS(1,:).^2+modGPS(2,:).^2),modGPS(3,:)]);
    maxdim=max(max(obsGPS')-min(obsGPS'));
    s=maxdim/max([maxobs,maxmodel])/2; % scale max vector length to half the max axis range
    
    figure('Position', [1, 1, 1200, 1000]);
     % Plot data displacements
     quiver(obsGPS(1,:), obsGPS(2,:), gps{i}.displacements(1,:)*s, gps{i}.displacements(2,:)*s, 0, ...
        'Color','k','LineWidth',1)
    hold on
    text(obsGPS(1,end), obsGPS(2,end)-2000,[num2str(scalebar*1000),' mm']) % Add scalebar
    quiver(obsGPS(1,:), obsGPS(2,:), zeros(size(gps{i}.displacements(3,:))), gps{i}.displacements(3,:)*s, 0, ...
        'Color','k','LineWidth',3,'MaxHeadSize',0.01,'Marker','s')

    % Plot modelled displacements
    quiver(obsGPS(1,1:end-1),obsGPS(2,1:end-1),modGPS(1,:)*s,modGPS(2,:)*s,0,'Color','r','LineWidth',1)
    quiver(obsGPS(1,1:end-1),obsGPS(2,1:end-1),zeros(size(modGPS(3,:))),modGPS(3,:)*s,0,'Color','r','LineWidth',1,'MaxHeadSize',0.01,'Marker','s')

    for i2 = 1:invpar.nModels % For each source model...
        if isempty(invpar.modelinsarID{i2}) || ismember(i, invpar.modelinsarID{i2}) % if model included for this GPS data set
            index1 = invResults.model.mIx(i2);
            switch invpar.model{i2}
                case {'MOGI','CDM','PCDM','FECM'}
                    plot(invResults.model.optimal(index1), invResults.model.optimal(index1+1),'b*')
                case {'DIKE','FAUL','DLOC','MDIK','MDLC'}
                    m=[invResults.optimalmodel{i2}];
                    drawmodel(m,'color','b','updipline','yes','projection','no')
            end
        end
    end



    axis equal; 
    ax = gca;
    grid on
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.LineWidth = 1.5;
    ax.GridLineStyle = '--';
    xlabel('X distance from local origin (m)')
    ylabel('Y distance from local origin (m)')
    title('GPS horizontal displacements (data:black - model:red)')
    %xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
    %ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);

    drawnow




    if saveflag=='y'
        print(gcf,[outputDir,'/Figures/GPS_Data_Model_',num2str(i)],'-dpng')
    end
end
%drawnow

