function plotInSARDataModelResidual(insar, geo, invpar, invResults, modelinput, saveName, fidHTML, saveflag, singlePlotFigures)

% Function to generate plot with comparison between InSAR data, model, and
% residuals
%
% Usage: plotInSARDataModelResidual(insar, geo, invpar, invResults, modelinput, saveName, fidHTML, saveflag)
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

% Create colormaps for plotting InSAR data
cmap.seismo = colormap_cpt('GMT_seis.cpt', 100); % GMT Seismo colormap for wrapped interferograms
cmap.redToBlue = colormap_cpt('polar.cpt', 100); % Red to Blue colormap for unwrapped interferograms

for i=1:length(insar)
    % Load and display DATA
    if isfield(insar{i},'plotPath')
        loadedData = load(insar{i}.plotPath); % load *.mat file
    else    
        loadedData = load(insar{i}.dataPath); % load *.mat file
    end
    
%     % Apply bounding box and remove data points outside it
%     iOutBox = find(loadedData.Lon<geo.boundingBox(1) | loadedData.Lon>geo.boundingBox(3) | loadedData.Lat>geo.boundingBox(2) | loadedData.Lat<geo.boundingBox(4));
%     if sum(iOutBox)>0
%         loadedData.Phase(iOutBox) = [];
%         loadedData.Lat(iOutBox) = [];
%         loadedData.Lon(iOutBox) = [];
%         loadedData.Heading(iOutBox) = [];
%         loadedData.Inc(iOutBox) = [];
%     end
%     
    convertedPhase = (loadedData.Phase/(4*pi))*insar{i}.wavelength;    % Convert phase from radians to m
    los = -convertedPhase;  % Convert phase from cm to Line-of-sigth displacement in m
    Heading = loadedData.Heading;
    Inc = loadedData.Inc;
    ll = [single(loadedData.Lon) single(loadedData.Lat)];   % Create Longitude and Latitude matrix
    xy = llh2local(ll', geo.referencePoint);    % Transform from geographic to local coordinates
    
    nPointsThis = size(ll,1);   % Calculate length of current data vector
    xy = double([[1:nPointsThis]',xy'*1000]);   % Add ID number column to xy matrix with local coordinates
    
    % Patch scattered data for faster plotting
    edge = round(min(abs(diff(xy(:,3)))))+2; % Size of patch set to minumum distance between points
    if edge < 50
        edge = 50;
    end
    xs = [xy(:,2)'; xy(:,2)'+edge; xy(:,2)'+edge; xy(:,2)']; % Coordinates of four vertex of patch
    ys = [xy(:,3)'; xy(:,3)'; xy(:,3)'+edge; xy(:,3)'+edge];
    
    % Extract filename to be included in figure name
    [pathstr,name,ext] = fileparts(insar{i}.dataPath);
    disp(name)
    % Display wrapped DATA interferogram at 5.6 cm wavelength
    figure('Position', [0.1, 0.1, 1000, 300]);
    if strcmpi(singlePlotFigures,'y') 
        ax1=gcf;
    else
        ax1 = subplot(1,3,1);
    end  
    plotInsarWrapped(xy,los, insar{i}.wavelength, cmap, 'DATA');
    colormap(ax1,cmap.redToBlue)
    % Add title
    title('data', 'Interpreter', 'none');


    % Create the 'GMT_input' directory if it doesn't exist
    gmtDir = fullfile(outputDir, 'GMT_input');
    if ~exist(gmtDir, 'dir')
        mkdir(gmtDir);
    end

    %merged txt file
    combined_data = [ll, los]; % This will create an n*3 matrix: [longitude, latitude, los]
    % Define the output file name
    outputFile = fullfile(gmtDir, [name '_raw_data.txt']);
    % Save the combined data into the specified file with tab delimiter
    writematrix(combined_data, outputFile, 'Delimiter', 'tab');

    % Calculate MODEL
    constOffset = 0;
    xRamp = 0;
    yRamp = 0;
    
    if i == 1
        if insar{i}.constOffset == 'y'
            constOffset = invResults.model.mIx(end);
            invResults.model.mIx(end) = invResults.model.mIx(end)+1;
        end
        if insar{i}.rampFlag == 'y'
            xRamp = invResults.model.mIx(end);
            yRamp = invResults.model.mIx(end)+1;
            invResults.model.mIx(end) = invResults.model.mIx(end)+2;
        end
    end
    
   if i > 1
        if insar{i}.constOffset == 'y'
            constOffset = invResults.model.mIx(end);
            invResults.model.mIx(end) = invResults.model.mIx(end)+1;
        end
        if insar{i}.rampFlag == 'y'
            xRamp = invResults.model.mIx(end);
            yRamp = invResults.model.mIx(end)+1;
            invResults.model.mIx(end) = invResults.model.mIx(end)+2;
        end
    end
       
    [modLos, UTot] = forwardInsarModel(insar{i},xy,invpar,invResults,modelinput,geo,Heading,Inc,constOffset,xRamp,yRamp,i); % Modeled InSAR displacements
 
    %save one time for ENU
    E = [ll, UTot(1,:)'];
    N = [ll, UTot(2,:)'];
    U = [ll, UTot(3,:)'];
    % This will create an n*3 matrix: [longitude, latitude, los]
    % Define the output file name
    savenameENU = erase(saveName, '.mat');
    outputFileE = fullfile(gmtDir, [savenameENU 'Emodel.txt']);
    outputFileN = fullfile(gmtDir, [savenameENU 'Nmodel.txt']);
    outputFileU = fullfile(gmtDir, [savenameENU 'Umodel.txt']);

        % Check if each file exists before saving, and write if not
    if ~exist(outputFileE, 'file')
        writematrix(E, outputFileE, 'Delimiter', 'tab');
    end
    
    if ~exist(outputFileN, 'file')
        writematrix(N, outputFileN, 'Delimiter', 'tab');
    end
    
    if ~exist(outputFileU, 'file')
        writematrix(U, outputFileU, 'Delimiter', 'tab');
    end



    % Display MODEL of data 
     if strcmpi(singlePlotFigures,'y') 
        ax2=figure('Position', [1, 1, 2400, 2000]);
    else
        ax2 = subplot(1,3,2);
    end   
    plotInsarWrapped(xy,modLos',insar{i}.wavelength,  cmap, 'MODEL');
    colormap(ax2,cmap.redToBlue)
    % Add title
    title('model_data', 'Interpreter', 'none');
        
    %merged txt file
    combined_data = [ll, modLos']; % This will create an n*3 matrix: [longitude, latitude, los]
    % Define the output file name
    outputFile = fullfile(gmtDir, [name '_model.txt']);
    % Save the combined data into the specified file with tab delimiter
    writematrix(combined_data, outputFile, 'Delimiter', 'tab');
    
%     % Display MODEL of offset
%     if strcmpi(singlePlotFigures,'y') 
%         ax3=figure('Position', [1, 1, 1200, 1000]);
%     else
%         ax3 = subplot(1,6,3);
%     end
%     plotInsarWrapped(xy,modoffset',insar{i}.wavelength, cmap, 'MODEL');
%     colormap(ax3,cmap.redToBlue)
%     caxis([-0.5 0.5])
%     % Add title
%     title('model_offset', 'Interpreter', 'none');
% 
% 
%     % Display MODEL of offset
%     if strcmpi(singlePlotFigures,'y') 
%         ax4=figure('Position', [1, 1, 1200, 1000]);
%     else
%         ax4 = subplot(1,6,4);
%     end
%     plotInsarWrapped(xy,modramp',insar{i}.wavelength, cmap, 'MODEL');
%     colormap(ax4,cmap.redToBlue)
%     caxis([-0.5 0.5])
%     % Add title
%     title('model_ramp', 'Interpreter', 'none');


    % Display RESIDUAL (data-model)
    residual = los-modLos';
    % Calculate the RMS of residuals
    rmsResidual = sqrt(nanmean(residual.^2));
    disp(['RMS Residual_model: ', num2str(rmsResidual)]);


    if strcmpi(singlePlotFigures,'y') 
        ax5=figure('Position', [1, 1, 1200, 1000]);
    else
        ax5 = subplot(1,3,3);
    end
    plotInsarWrapped(xy,residual, insar{i}.wavelength, cmap, 'RESIDUAL');
    caxis([-1 1])
    colormap(ax5,cmap.redToBlue)
    % Add title
    title('Res(data - m_data)', 'Interpreter', 'none');



%     % Display RESIDUAL (data-model-offset-ramp)
%     residual = los-modLos'-modoffset'-modramp';
%     % Calculate the RMS of residuals
%     rmsResidual = sqrt(mean(residual.^2));
%     disp(['RMS Residual_off: ', num2str(rmsResidual)]);
% 
% 
%     if strcmpi(singlePlotFigures,'y') 
%         ax6=figure('Position', [1, 1, 1200, 1000]);
%     else
%         ax6 = subplot(1,6,6);
%     end
%     plotInsarWrapped(xy,residual, insar{i}.wavelength, cmap, 'RESIDUAL');
%     caxis([-1 1])
%     colormap(ax6,cmap.redToBlue)
%     % Add title
%     title('Res(data - m_data - m_offset - m_ramp)', 'Interpreter', 'none');
%  


    if saveflag=='y'
        if strcmpi(singlePlotFigures,'y') 
            img = getframe(ax1);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Wrap_Data_',name,'.png']);
            img = getframe(ax2);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Unwrap_Data_',name,'.png']);            
            img = getframe(ax3);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Wrap_Model_',name,'.png']);            
            img = getframe(ax4);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Unwrap_Model_',name,'.png']);            
            img = getframe(ax5);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Wrap_Residual_',name,'.png']);            
            img = getframe(ax6);
            imwrite(img.cdata,[outputDir,'/Figures/InSAR_Unwrap_Residual_',name,'.png']);            
        else
            %img = getframe(gcf);
            %imwrite(img.cdata,[outputDir,'/Figures/InSAR_Data_Model_Residual_',name,'.png']);
            print('-dpng', [outputDir,'/Figures/InSAR_Data_Model_Residual_',name,'.png'])
        end        
        
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Comparison InSAR Data - Model - Residual</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/InSAR_Data_Model_Residual_',name,'.png','" alt="HTML5 Icon">']);
    end
    
    insar_ix = i;
end

end % end function

function []=plotSourceOutlines(invpar,invResults)
    
hold on
for i2 = 1:invpar.nModels % For each source model...
  if isempty(invpar.modelinsarID{i2}) || ismember(i, invpar.modelinsarID{i2}) % if model included for this InSAR data set      index1 = invResults.model.mIx(i2);
      index1 = invResults.model.mIx(i2);
      switch invpar.model{i2}
        case {'MOGI','CDM','PCDM','FECM'}
            plot(invResults.model.optimal(index1)/1000, invResults.model.optimal(index1+1)/1000,'k*','markersize',14)
        case {'DIKE','FAUL','DLOC','MDIK','MDLC'}
            m=[invResults.optimalmodel{i2}];
            m([1:3,6,7],:)=m([1:3,6,7],:)/1000;
            drawmodel(m,'color','k','updipline','yes','projection','no')     
        case {'SILL'}
            m=[invResults.optimalmodel{i2}];
            m([1:3,6,7],:)=m([1:3,6,7],:)/1000;
            drawmodel(m,'color','k','updipline','no','projection','no')     
      end
  end
end
hold off

end % end function