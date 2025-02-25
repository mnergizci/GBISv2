function [ULos_model, UTot] = forwardInSARModel(insar,xy,invpar,invResults,modelInput,geo,Heading,Inc,constOffset,xRamp,yRamp,insarID)

% Function to generate forward model for InSAR displacements using optimal
% source parameters
%
% Usage: ULos = forwardInSARModel(insar,xy,invpar,invResults,modelInput,geo,Heading,Inc,constOffset,xRamp,yRamp)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
% =========================================================================
%%

nu = modelInput.nu;
xy = [xy(:,2:3)'; zeros(1,length(xy(:,1)))];
UTot = zeros(3,length(xy(1,:)));   % Initialise matrix of modeled displacements (3 x number of observation points)

for i = 1:invpar.nModels % For each source model...
    index1 = invResults.model.mIx(i);
    if isempty(invpar.modelinsarID{i}) || ismember(insarID, invpar.modelinsarID{i}) % if model included for this InSAR data set
        switch invpar.model{i}
            case 'MOGI'
                mFunc{i} = invResults.model.optimal(index1:index1+3); % Parameters to invert for.
                U = mogi(mFunc{i},xy,nu);
            case 'MCTG'
                mFunc{i} = invResults.model.optimal(index1:index1+4);
                U = mctigueSource(mFunc{i},xy(1:2,:),nu);
            case 'YANG'
                mFunc{i} = invResults.model.optimal(index1:index1+7);
                U = yangSource(mFunc{i},xy,nu);
            case 'PENN'
                mFunc{i}= invResults.model.optimal(index1:index1+4);
                U = pennySource(mFunc{i},xy,nu);
            case 'SILL'
                mFunc{i}=[invResults.model.optimal(index1:index1+2);0; ...
                    invResults.model.optimal(index1+3:index1+5);0;0; ...
                    invResults.model.optimal(index1+6)]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'DIKE'
                mFunc{i}=[invResults.model.optimal(index1:index1+6);0;0;...
                    invResults.model.optimal(index1+7)]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'FAUL'
                mFunc{i}=[invResults.model.optimal(index1:index1+8);0]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'DLOC'
                mFunc{i}=[invResults.model.optimal(index1:index1+9)]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'MDIK'
                mFunc{i}=[invResults.optimalmodel{i}]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'MDLC'
                mFunc{i}=[invResults.optimalmodel{i}]; 
                U = disloc(mFunc{i},xy(1:2,:),nu); 
            case 'HING'
                mFunc{i}=invResults.model.optimal(index1:index1+10);
                U = hingedDikes(mFunc{i},xy(1:2,:),nu);
            case 'PCDM'
                mFunc{i}=invResults.model.optimal(index1:index1+8); 
                [ue,un,uv] = pCDM(xy(1,:),xy(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:),mFunc{i}(4,:),mFunc{i}(5,:),mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),nu);
                U=[ue';un';uv'];
            case 'CDM'
                mFunc{i}=invResults.model.optimal(index1:index1+9);
                [ue,un,uv] = CDM(xy(1,:),xy(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:),mFunc{i}(4,:),mFunc{i}(5,:), ...
                    mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),mFunc{i}(10,:),nu);
                U=[ue';un';uv'];
            case 'FECM'
                mu = modelInput.mu; % shear modulus
                lambda=2*mu*nu/(1-2*nu); % Lamé first constant
                mFunc{i}=invResults.model.optimal(index1:index1+9);
                [ue,un,uv,dV,DV] = fECM(xy(1,:),xy(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:),mFunc{i}(4,:),mFunc{i}(5,:), ...
                   mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),mFunc{i}(10,:),mu,lambda,'C');
                U=[ue';un';uv'];
        end
        UTot = UTot + U; % Calculate total displacement from sum of displacement from each source
    end
end

ULos_model = [];
% ULos_offset = [];
% ULos_ramp = []; % Initialize ramp to 0

UEast = -cosd(Heading).* sind(Inc); % East unit vector
UNorth = sind(Heading).* sind(Inc); % North unit vector
UVert = cosd(Inc); % Vertical unit vector

ULos_model = UEast'.* UTot(1,:) + ...
    UNorth'.* UTot(2,:) + ... % Convert to line of sight displacement
    UVert'.* UTot(3,:);
% 
% if insar.constOffset == 'y'
%     ULos_offset = invResults.model.optimal(constOffset);  % Add constant offset
%     disp(ULos_offset)
%     ULos_offset = ULos_offset * ones(size(ULos_model)); % Ensure correct size
% end
% 
% if insar.constOffset == 'n'
%     ULos_offset = 0;  % Add constant offset
%     disp(ULos_offset)
%     ULos_offset = ULos_offset * ones(size(ULos_model)); % Ensure correct size
% end
% 
% 
% if insar.rampFlag == 'y'
%     ULos_ramp = invResults.model.optimal(xRamp) * xy(1,:) + ...
%                 invResults.model.optimal(yRamp) * xy(2,:); % Add ramp
% end
% 
% if insar.rampFlag == 'n'
%     ULos_ramp = 0;
%     ULos_ramp = ULos_ramp * ones(size(ULos_model)); % Ensure correct size
% end

end
