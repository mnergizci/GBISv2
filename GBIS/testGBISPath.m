function testGBISPath()

% Test if paths to GBIS are set correctly
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

disp(' ')
which GBISrun
which loadInsarData
which mogi
which PlotInsarWrapped
which generateFinalReport
which runInversion
which local2llh
which fitVariogram
disp(' ')
disp('If you see this message, you are ready to go!')