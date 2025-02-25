function [] = GBISclean(inputDir)

% Delete output directory
%
% Usage:  GBISclean('directoryName')
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

%% Check for correct input
if nargin == 0
    error('Directory to remove not specified.')
end

%% Display dialog window to confirm deletion
choice = questdlg(['Do really you want to delete the directory "',inputDir,'"?'], 'Warning!', 'Yes', 'No','Yes');

switch choice
    case 'Yes'
        rmdir(inputDir,'s')
    case 'No'
        return
end