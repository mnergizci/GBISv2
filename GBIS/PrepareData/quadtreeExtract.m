function [values] = quadtreeExtract(xy,val,xlims,ylims)

% Function to extract values based on quadtree results. Modified after Gonzalez (2015).
%
% Usage: [values]=quadtreeExtract(xy,val,xlims,ylims)
% Input Parameters:
%   xy    : Mx3 array of pixels positions
%   val   : Mx1 array of the value of the pixels
%   xlims : nx2 array of x-coordinate limits for each selected polygon
%   ylims : nx2 array of y-coordinate limits for each selected polygon
%
% Output Parameters:
%   values: n   mean value inside each polygon
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

values    = zeros(1,length(xlims));
  
  % Loop over all quadtree subdivisions
  for i=1:length(xlims)
    in = find( xy(:,2)<=xlims(i,3) & xy(:,2)>=xlims(i,1) & xy(:,3)<=ylims(i,2) & xy(:,3)>=ylims(i,1));
    values(:,i) = mean(val(in));
  end

end