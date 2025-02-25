function generateFinalReport(invResFile,burning,singlePlotFigures)

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


%% Print results to text file
% Print optimal, mean, median, 2.5th and 97.5th percentiles 
format shortG

txtName = ['/summary',saveName(7:end-4),'.txt']; % name of text file
fileID = fopen([outputDir,txtName],'w');
fprintf(fileID,'%s\r\n','GBIS');
fprintf(fileID,'%s\r\n',['Summary for ',saveName(1:end-4)]);
fprintf(fileID,'%s\r\n',['Number of iterations: ',num2str(invpar.nRuns)]);
fprintf(fileID,'%s\r\n',['Burn-in (n. of iterations): ',num2str(burning)]);
fprintf(fileID,'%s\r\n','=================================================================================================');
nParamDist=0; % number of parms for distributed slip/opening
if isfield(invResults.model,'modelName')
    fprintf(fileID,'%s\r\n','MODEL     PARAM.          OPTIMAL            MEAN          Median            2.5%           97.5%');
    for i = 1:nParam-1
        switch invResults.model.modelName{i}(1:end-1)
            case {'MDisloc', 'MDike'}
                nParamDist=nParamDist+1;
            otherwise
                fprintf(fileID,'%s\t %7s\t %8g\t %8g\t %8g\t %8g\t %8g\r\n', ...
                char(invResults.model.modelName(i)),char(invResults.model.parName(i)), round([invResults.model.optimal(i), ...
                mean(invResults.mKeep(i, burning:end-blankCells)), median(invResults.mKeep(i, burning:end-blankCells)), ...
                prctile(invResults.mKeep(i, burning:end-blankCells),2.5), prctile(invResults.mKeep(i, burning:end-blankCells),97.5)],4,'significant'));
                if strcmp(invResults.model.parName(i),'Opening')
                    fprintf(fileID,'%s\t %7s\t %8g\t %8g\t %8g\t %8g\t %8g\r\n', ...
                    char(invResults.model.modelName(i)),'DV', round([invResults.model.optimal(i), ...
                    mean(invResults.mKeep(i, burning:end-blankCells)), median(invResults.mKeep(i, burning:end-blankCells)), ...
                    prctile(invResults.mKeep(i, burning:end-blankCells),2.5), prctile(invResults.mKeep(i, burning:end-blankCells),97.5)],4,'significant'));            
                end
        end
    end
else % for backwards compatibility
    fprintf(fileID,'%s\r\n','PARAMETER          OPTIMAL           MEAN          Median            2.5%           97.5%');
    for i = 1:nParam-1
        fprintf(fileID,'%12s\t %8g\t %8g\t %8g\t %8g\t %8g\r\n', ...
        char(invResults.model.parName(i)), round([invResults.model.optimal(i), ...
        mean(invResults.mKeep(i, burning:end-blankCells)), median(invResults.mKeep(i, burning:end-blankCells)), ...
        prctile(invResults.mKeep(i, burning:end-blankCells),2.5), prctile(invResults.mKeep(i, burning:end-blankCells),97.5)],4,'significant'));
    end
end

% Display text on screen
clc
type([outputDir,txtName])

%% Create report html file
htmlName = ['/report',saveName(7:end-4),'.html']; % HTML file name

% Print header and summary to file
fidHTML = fopen([outputDir,htmlName],'w');
fprintf(fidHTML, '%s\r\n', '<!DOCTYPE html>');
fprintf(fidHTML, '%s\r\n', '<html>');
fprintf(fidHTML, '%s\r\n', '<head>');
fprintf(fidHTML, '%s\r\n', ['<H1>GBIS Final report for <i>', saveName(1:end-4),'</i></H1>']);
fprintf(fidHTML, '%s\r\n', ['<H3>Results file: <i>',outputDir,'/',saveName,'</i></H3>']);
fprintf(fidHTML, '%s\r\n', ['<p>Number of iterations: ', num2str(invpar.nRuns),'</p>']);
fprintf(fidHTML, '%s\r\n', ['<p>Burning time (n. of iterations from start): ', num2str(burning),'</p>']);
fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Model parameters</H3>');
fprintf(fidHTML, '%s\r\n', '<style>');
fprintf(fidHTML, '%s\r\n', 'table {font-family: arial, sans-serif; border-collapse: collapse; width:100%%;}');
fprintf(fidHTML, '%s\r\n', 'td, th {border: 1px solid #dddddd;text-align: right;padding: 8px;}');
fprintf(fidHTML, '%s\r\n', 'tr:nth-child(even) {background-color: #bbb;}');
fprintf(fidHTML, '%s\r\n', '</style>');
fprintf(fidHTML, '%s\r\n', '</head>');
fprintf(fidHTML, '%s\r\n', '<body>');
fprintf(fidHTML, '%s\r\n', '<table>');
if isfield(invResults.model,'modelName')
    fprintf(fidHTML, '%s\r\n', '<tr> <th>Model</th> <th>Parameter</th> <th>Optimal</th> <th>Mean</th> <th>Median</th> <th>2.5%</th> <th>97.5%</th></tr>');
    for i = 1:nParam-1
        fprintf(fidHTML, '%s\r\n',['<tr> <td>',char(invResults.model.modelName(i)),...
        '</td> <td>',char(invResults.model.parName(i)),...
        '</td> <td>', num2str(invResults.model.optimal(i), '%.2f'), ...
        '</td> <td>', num2str(mean(invResults.mKeep(i,burning:end-blankCells)), '%.2f'), ...
        '</td> <td>', num2str(median(invResults.mKeep(i,burning:end-blankCells)), '%.2f'),...
        '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),2.5), '%.2f'), ...
        '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),97.5), '%.2f'), ...
        '</td> </tr>']);    
    end
else % for backwards compatibility
    fprintf(fidHTML, '%s\r\n', '<tr> <th>Parameter</th> <th>Optimal</th> <th>Mean</th> <th>Median</th> <th>2.5%</th> <th>97.5%</th></tr>');
    for i = 1:nParam-1
        fprintf(fidHTML, '%s\r\n',['</td> <td>',char(invResults.model.parName(i)),...
        '</td> <td>', num2str(invResults.model.optimal(i), '%.2f'), ...
        '</td> <td>', num2str(mean(invResults.mKeep(i,burning:end-blankCells)), '%.2f'), ...
        '</td> <td>', num2str(median(invResults.mKeep(i,burning:end-blankCells)), '%.2f'),...
        '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),2.5), '%.2f'), ...
        '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),97.5), '%.2f'), ...
        '</td> </tr>']);
    end
end
fprintf(fidHTML, '%s\r\n', '</table> </body>');



%% Plot comparison betweem data, model, and residual

% Optional
%choice = questdlg('Do you want to compare DATA MODEL and RESIDUAL?', 'Plot?', 'Yes', 'No','Yes');
choice = 'No';
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
end

%% Plot distributed opening if there's a distributed dike model



%% Plot convergence of all parameters
figure('Position', [1, 1, 1200, 1000]);
plot(1:100:length(invResults.PKeep(1,:))-blankCells, invResults.PKeep(1,1:100:end-blankCells),'r.','MarkerSize',10) % Plot one point every 100 iterations
title('Probability Exponent')

img = getframe(gcf);
imwrite(img.cdata, [outputDir,'/Figures/Probability.png']);

fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Proability exponent plot</H3>');
fprintf(fidHTML, '%s\r\n', ['<img src="Figures/Probability.png','" alt="HTML5 Icon">']);
        
%choice = questdlg('Do you want to plot convergence figures?', 'Plot?', 'Yes', 'No','Yes');
choice = 'Yes';
switch choice
    case 'Yes'
        iPlot = 0; % Counter for subplots in the single image
        separatePlotIndex = 1; % Counter for separate saved images
        
        % Create a figure for combined plots
        figure('Position', [1, 1, 1200, 1000]);
        % Number of parameters to select
        numRandomPlots = min(50, nParam-1); % Select 50 or fewer parameters if fewer exist
        randomIndices = randperm(nParam-1, numRandomPlots); % Randomly select 50 unique indices

        for idx = randomIndices % Loop only through the randomly selected indices
            i = idx; % Use the randomly selected index
            switch invResults.model.modelName{i}(1:end-1)
                case {'MDisloc', 'MDike'}
                    % Save distributed parameters as separate figures
                    separateFigure = figure('Position', [100, 100, 800, 600]);
                    plot(1:100:length(invResults.mKeep(1,:))-blankCells, invResults.mKeep(i,1:100:end-blankCells), 'r.');
                    title([invResults.model.modelName{i}, ' ', invResults.model.parName{i}]);
                    xlabel('Iterations');
                    ylabel('Parameter Value');
                    
                    % Save the figure
                    saveFileName = [outputDir, '/Figures/converge', num2str(separatePlotIndex), '.png'];
                    saveas(separateFigure, saveFileName);
                    close(separateFigure); % Close the figure
                    
                    % Increment the counter for separate plots
                    separatePlotIndex = separatePlotIndex + 1;

                    % Add the plot to the HTML report
                    fprintf(fidHTML, '%s\r\n', '<hr>');
                    fprintf(fidHTML, '%s\r\n', ['<H3>Convergence plot for ', invResults.model.modelName{i}, ' ', invResults.model.parName{i}, '</H3>']);
                    fprintf(fidHTML, '%s\r\n', ['<img src="Figures/converge', num2str(separatePlotIndex - 1), '.png" alt="Convergence Plot">']);
%                 otherwise
%                     % Add non-distributed parameters to the single image
%                     iPlot = iPlot + 1;
%                     subplot(round((numRandomPlots-nParamDist)/3), 3, iPlot);
%                     plot(1:100:length(invResults.mKeep(1,:))-blankCells, invResults.mKeep(i,1:100:end-blankCells), 'r.');
%                     title([invResults.model.modelName{i}, ' ', invResults.model.parName{i}]);
            end
        end

        % Save the combined figure
        combinedSaveFileName = [outputDir, '/Figures/ConvergenceCombined.png'];
        saveas(gcf, combinedSaveFileName);
        close(gcf); % Close the combined figure
        
        % Add combined figure to HTML report
        fprintf(fidHTML, '%s\r\n', '<hr>');
        fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots for non-distributed parameters</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/ConvergenceCombined.png" alt="Combined Convergence Plot">']);


%% Plot histograms and optimal values
%choice = questdlg('Do you want to plot the individual PDFs?', 'Plot?', 'Yes', 'No','Yes');
choice = 'Yes';
switch choice
    case 'Yes'
        separatePlotIndex = 1; % Counter for separate plots
        
        % Determine random parameters for histogram plotting
        randomIndices = randperm(nParam-1, min(50, nParam-1)); % Randomly select up to 50 parameters
        
        % Create a figure for combined histograms
        combinedFigure = figure('Position', [1, 1, 1200, 1000]);
        iPlot = 0; % Counter for subplot index
        
        for i = randomIndices
            switch invResults.model.modelName{i}(1:end-1)
                case {'MDisloc', 'MDike'}
                    % Create separate histogram for distributed parameters
                    separateFigure = figure('Position', [100, 100, 800, 600]);
                    xMin = mean(invResults.mKeep(i, burning:end-blankCells)) - 4 * std(invResults.mKeep(i, burning:end-blankCells));
                    xMax = mean(invResults.mKeep(i, burning:end-blankCells)) + 4 * std(invResults.mKeep(i, burning:end-blankCells));
                    bins = xMin: (xMax-xMin)/50: xMax;
                    h = histogram(invResults.mKeep(i, burning:end-blankCells), bins, 'EdgeColor', 'none', 'Normalization', 'count');
                    hold on;
                    topLim = max(h.Values);
                    plot([invResults.model.optimal(i), invResults.model.optimal(i)], [0, topLim + 10000], 'r-'); % Plot optimal value
                    ylim([0 topLim + 10000]);
                    title([invResults.model.modelName{i}, ' ', invResults.model.parName{i}]);
                    xlabel('Parameter Value');
                    ylabel('Frequency');
                    
                    % Save the separate figure
                    saveFileName = [outputDir, '/Figures/histogram', num2str(separatePlotIndex), '.png'];
                    saveas(separateFigure, saveFileName);
                    close(separateFigure); % Close the figure
                    
                    % Increment the counter for separate plots
                    separatePlotIndex = separatePlotIndex + 1;

%                 otherwise
%                     % Combine non-distributed parameter histograms in one figure
%                     iPlot = iPlot + 1;
%                     subplot(round(min(50, nParam-1)/3), 3, iPlot); % Determine subplot position
%                     xMin = mean(invResults.mKeep(i, burning:end-blankCells)) - 4 * std(invResults.mKeep(i, burning:end-blankCells));
%                     xMax = mean(invResults.mKeep(i, burning:end-blankCells)) + 4 * std(invResults.mKeep(i, burning:end-blankCells));
%                     bins = xMin: (xMax-xMin)/50: xMax;
%                     h = histogram(invResults.mKeep(i, burning:end-blankCells), bins, 'EdgeColor', 'none', 'Normalization', 'count');
%                     hold on;
%                     topLim = max(h.Values);
%                     plot([invResults.model.optimal(i), invResults.model.optimal(i)], [0, topLim + 10000], 'r-'); % Plot optimal value
%                     ylim([0 topLim + 10000]);
%                     title(invResults.model.parName(i));
            end
        end
        
        % Save the combined histograms
        combinedSaveFileName = [outputDir, '/Figures/PDFs.png'];
        saveas(combinedFigure, combinedSaveFileName);
        close(combinedFigure); % Close the combined figure
        
        % Add combined histograms and separate plots to HTML report
        fprintf(fidHTML, '%s\r\n', '<hr>');
        fprintf(fidHTML, '%s\r\n', '<H3>Combined Histograms for Non-Distributed Parameters</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/PDFs.png" alt="Combined Histograms">']);
        
        % Add separate histograms to HTML report
        for idx = 1:separatePlotIndex-1
            fprintf(fidHTML, '%s\r\n', '<hr>');
            fprintf(fidHTML, '%s\r\n', ['<H3>Histogram for Distributed Parameter ', num2str(idx), '</H3>']);
            fprintf(fidHTML, '%s\r\n', ['<img src="Figures/histogram', num2str(idx), '.png" alt="Histogram">']);
        end


%% Plot joint probabilities
%choice = questdlg('Do you want to plot the joint PDFs?', 'Plot?', 'Yes', 'No','Yes');
choice = 'False';
switch choice
    case 'Yes'
        
        figure('Position', [1, 1, 1200, 1000]);
        plotmatrix_lower(invResults.mKeep(1:nParam-1,burning:end-blankCells)','contour');
        img1 = getframe(gcf);
        imwrite(img1.cdata,[outputDir,'/Figures/JointProbabilities.png']);
        
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<hr>');
        fprintf(fidHTML, '%s\r\n', '<H3>Joint probabilities</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/JointProbabilities.png','" alt="HTML5 Icon">']);
        
    case 'No'
end

    case 'No'
end

    case 'No'
end
fclose(fidHTML);






