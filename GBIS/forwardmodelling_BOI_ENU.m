function forwardmodelling_BOI(invResFile,burning,singlePlotFigures)

% Function to generate summary text file, results plots, and html report
%
% Usage: generateFinalReport(invResFile,burn_in)
% Input parameters:
%           invResFile: path and name to file with final results of
%                       inversion 
%                       (e.g.,'VolcanoExercise/invert_1__MOGI_DIKE.mat')
%           burn_in:    number of iterations to ignore in pdf histogram plot
%                       and in computation of mean/median/confidence interval
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
if nargin < 2
    disp('!!!!!!! Not enough input parameters. !!!!!!!!!')
    return;
end

if nargin < 3
    singlePlotFigures='n';
end

if ~strcmpi(invResFile(end-3:end),'.mat')
    dirfile=dir([invResFile,'/invert*.mat']);
    invResFile=[invResFile,'/',dirfile.name];
else
    dirfile=dir(invResFile);
end
    

clear  outputDir  % Clear global variables
global outputDir  % Set global variables
disp(outputDir)
outputDir=dirfile.folder;

warning('off','all')

load(invResFile);   % load results

%%% here just to get the full BOI from different models! 
% Get the base path (without the last directory)
basePath1 = fileparts(invResFile);
basePath = fileparts(fileparts(basePath1));
% Define the name of the folder and the .mat file
folderName = 'invert_1_2_3_4_5_6_7_8_9_Z';
matFileName = 'invert_1_2_3_4_5_6_7_8_9_Z.mat';

% Create the full data path by merging basePath, folderName, and matFileName
fulldatapath = fullfile(basePath, folderName, matFileName);

% Load the 'insar' variable from the .mat file
load(fulldatapath, 'insar');


%%redefine the full path of the insar to make forward modelling
for i = 1:length(insar)
    % Get the original path
    originalPath = insar{i}.dataPath;
    
    % Extract the base path (up to '.../Data')
    basePath = fileparts(fileparts(originalPath)); % This goes two levels up from the file path
    
    % Define the new folder name
    newFolder = 'TR_full';
    
    % Extract the base name of the file (excluding '_rr.mat')
    [~, fileBase, ~] = fileparts(originalPath);
    fileBase = erase(fileBase, '_rr'); % Remove '_rr' if it exists
    
    % Create the new file path
    newFilePath = fullfile(basePath, newFolder, [fileBase '.mat']);
    
    % Display the new full path
    disp(newFilePath);
    
    % (Optional) Store the new path back into the structure if needed
    insar{i}.plotPath = newFilePath;
end


% Create colormaps for plotting InSAR data
cmap.Seismo = colormap_cpt('GMT_seis.cpt', 100); % GMT Seismo colormap for wrapped interferograms
cmap.RnB = colormap_cpt('polar.cpt', 100); % Red to Blue colormap for unwrapped interferograms

nParam = length(invResults.mKeep(:,1)); % Number of model parameters


% Number of empty cells at the end of mKeep and pKeep
if invpar.nRuns < 10000
    blankCells = 999; 
else
    blankCells = 9999;
end

htmlName = ['/report',saveName(7:end-4),'.html']; % HTML file name

% Print header and summary to file
fidHTML = fopen([outputDir,htmlName],'w');


%% Plot comparison betweem data, model, and residual

% Optional
%choice = questdlg('Do you want to compare DATA MODEL and RESIDUAL?', 'Plot?', 'Yes', 'No','Yes');
choice = 'Yes';
switch choice
    case 'Yes'
        % Plot GPS data, model
        
        if exist('gps','var')
            plotGPSDataModel(gps,geo,invpar, invResults, modelInput, saveName, 'y')
                        
            % Add image to html report
            fprintf(fidHTML, '%s\r\n', '<hr>');
            fprintf(fidHTML, '%s\r\n', '<H3>Comparison GPS Data vs. Model</H3>');
            fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model.png','" alt="HTML5 Icon">']);
        end
        
        % Plot InSAR data, model, residual
        if exist('insarDataCode')
            plotInSARDataModelResidual(insar, geo, invpar, invResults, modelInput, saveName, fidHTML, 'y',singlePlotFigures)
        end
    case 'No'
        return
end

% %% Plot distributed opening if there's a distributed dike model
% 
% 
% 
% %% Plot convergence of all parameters
% figure('Position', [1, 1, 1200, 1000]);
% plot(1:100:length(invResults.PKeep(1,:))-blankCells, invResults.PKeep(1,1:100:end-blankCells),'r.','MarkerSize',10) % Plot one point every 100 iterations
% title('Probability Exponent')
% 
% img = getframe(gcf);
% imwrite(img.cdata, [outputDir,'/Figures/Probability.png']);
% 
% fprintf(fidHTML, '%s\r\n', '<hr>');
% fprintf(fidHTML, '%s\r\n', '<H3>Proability exponent plot</H3>');
% fprintf(fidHTML, '%s\r\n', ['<img src="Figures/Probability.png','" alt="HTML5 Icon">']);
%         
% %choice = questdlg('Do you want to plot convergence figures?', 'Plot?', 'Yes', 'No','Yes');
% choice = 'False';
% switch choice
%     case 'Yes'
% 
%         figure('Position', [1, 1, 1200, 1000]);
%         iPlot=0;
%         for i = 1:nParam-1
%             switch invResults.model.modelName{i}(1:end-1)
%                 case {'MDisloc', 'MDike'}
%                 otherwise
%                     iPlot=iPlot+1;
%                     subplot(round((nParam-nParamDist)/3),3,iPlot)    % Determine poistion in subplot
%                     plot(1:100:length(invResults.mKeep(1,:))-blankCells, invResults.mKeep(i,1:100:end-blankCells),'r.') % Plot one point every 100 iterations
%                     title([invResults.model.modelName(i),' ',invResults.model.parName(i)])
%             end
%         end
%         % Save image as png
%         img = getframe(gcf);
%         imwrite(img.cdata, [outputDir,'/Figures/Convergence.png']);
%         
%         % Add image to html report
%         fprintf(fidHTML, '%s\r\n', '<hr>');
%         fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots</H3>');
%         fprintf(fidHTML, '%s\r\n', ['<img src="Figures/Convergence.png','" alt="HTML5 Icon">']);
% 
% 
% 
% %% Plot histograms and optimal values
% %choice = questdlg('Do you want to plot the individual PDFs?', 'Plot?', 'Yes', 'No','Yes');
% choice = 'False';
% switch choice
%     case 'Yes'
% 
%         figure('Position', [1, 1, 1200, 1000]);
%         for i = 1:nParam-1
%             subplot(round(nParam/3),3,i) % Determine poistion in subplot
%             xMin = mean(invResults.mKeep(i,burning:end-blankCells))-4*std(invResults.mKeep(i,burning:end-blankCells));
%             xMax = mean(invResults.mKeep(i,burning:end-blankCells))+4*std(invResults.mKeep(i,burning:end-blankCells));
%             bins = xMin: (xMax-xMin)/50: xMax;
%             h = histogram(invResults.mKeep(i,burning:end-blankCells),bins,'EdgeColor','none','Normalization','count');
%             hold on
%             topLim = max(h.Values);
%             plot([invResults.model.optimal(i),invResults.model.optimal(i)],[0,topLim+10000],'r-') % Plot optimal value
%             ylim([0 topLim+10000])
%             title(invResults.model.parName(i))
%         end
%         % Save image as png
%         img = getframe(gcf);
%         imwrite(img.cdata,[outputDir,'/Figures/PDFs.png']);
%         
%         % Add image to html report
%         fprintf(fidHTML, '%s\r\n', '<hr>');
%         fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Model parameters posterior probabilities and optimal values</H3>');
%         fprintf(fidHTML, '%s\r\n', ['<img src="Figures/PDFs.png','" alt="HTML5 Icon">']);
% 
% 
% 
% %% Plot joint probabilities
% %choice = questdlg('Do you want to plot the joint PDFs?', 'Plot?', 'Yes', 'No','Yes');
% choice = 'False';
% switch choice
%     case 'Yes'
%         
%         figure('Position', [1, 1, 1200, 1000]);
%         plotmatrix_lower(invResults.mKeep(1:nParam-1,burning:end-blankCells)','contour');
%         img1 = getframe(gcf);
%         imwrite(img1.cdata,[outputDir,'/Figures/JointProbabilities.png']);
%         
%         % Add image to html report
%         fprintf(fidHTML, '%s\r\n', '<hr>');
%         fprintf(fidHTML, '%s\r\n', '<H3>Joint probabilities</H3>');
%         fprintf(fidHTML, '%s\r\n', ['<img src="Figures/JointProbabilities.png','" alt="HTML5 Icon">']);
%         
%     case 'No'
% end
% 
%     case 'No'
% end
% 
%     case 'No'
% end
% fclose(fidHTML);
% 





