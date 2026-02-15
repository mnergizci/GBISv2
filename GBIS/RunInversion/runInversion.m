function results = runInversion(geo, gps, enu, insar, invpar, model, modelInput, obs, nObs)

% Function that runs the MCMC Bayesian inversion
%
% Usage: results = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs)
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

global outputDir  % Set global variables

%% Set starting model parameters
nu = modelInput.nu; % Poisson's ratio

model.trial = model.m; % First trial using starting parameters

model.range = model.upper - model.lower; % Calculate range of m parameters
% Check that starting model is within bounds
if sum(model.m > model.upper)>0 || sum(model.m < model.lower)>0
    disp('Parameter#  Lower_Bound  Sart_Model  Upper_Bound')
    ix = find(model.m > model.upper | model.m < model.lower);
    for i=1:length(ix)
        fprintf('    %d      %f %f %f\n', ix(i), model.lower(ix(i)), model.m(ix(i)), model.upper(ix(i)))
    end
    error('Starting model is out of bounds')
end

nModel = length(model.m);    % Number of model parameters to invert for
probTarget = 0.5^(1/nModel); % Target probability used in sensitivity test
probSens = zeros(nModel,1);  % Initialise vector for sensitivity test

iKeep = 0; % Initialiase kept iterations counter
iReject = 0;    % Initialise rejected iterations counter
iKeepSave = iKeep;  % Initialise saving schedule for kept iterations
iRejectSave = iReject; % Initialise saving schedule for rejected iterations

mKeep= zeros(nModel,invpar.nSave,'single');   % Initialise matrix of model parameters to keep
PKeep= zeros(1,invpar.nSave,'single');   % Initialise probability vector

POpt = -1e99; % Set initial optimal probability
T = invpar.TSchedule(1); % Set first temperature
iTemp = 0;  % Initialise temperature schedule index (for initial Simulated Annealing)
nTemp = length(invpar.TSchedule); % Number of temperature steps (for initial Simulated Annealing)

i_mdike=0;
i_mdisloc=0;
for i = 1:invpar.nModels % For each source model...
    index1 = model.mIx(i);
    switch invpar.model{i}
        case 'MDIK'
            i_mdike=i_mdike+1;
            mdike=load(modelInput.mdike{i_mdike}.patchfile);
            mdike.m([8,9],:)=0; % Set slip to 0
            mdike.m(10,:)=1; % Set openning to 1
            mFunc{i}=mdike.m;
            for i2 = 1:size(mdike.m,2)
                Uunit{i}{i2} = disloc(mdike.m(:,i2),obs(1:2,:),nu); % displacements for unit slip
            end
        case 'MDLC'
            i_mdisloc=i_mdisloc+1;
            mdisloc=load(modelInput.mdisloc{i_mdisloc}.patchfile);
            mdisloc.m([8:10],:)=0; % Set slip and opening to 0
            mFunc{i}=mdisloc.m;
            ssdsop{i}=modelInput.mdisloc{i_mdisloc}.ss_ds_op;
            if ssdsop{i}(1)
                mdisloc.m(8,:)=1; % Set ss to 1
                for i2 = 1:size(mdisloc.m,2)
                    Uunit{i}{1}{i2} = disloc(mdisloc.m(:,i2),obs(1:2,:),nu); % displacements for unit slip
                end
                mdisloc.m(8,:)=0; % Set ss to 0
            end
            if ssdsop{i}(2)
                mdisloc.m(9,:)=1; % Set ds to 1
                for i2 = 1:size(mdisloc.m,2)
                    Uunit{i}{2}{i2} = disloc(mdisloc.m(:,i2),obs(1:2,:),nu); % displacements for unit slip
                end
                mdisloc.m(9,:)=0; % Set ds to 0
            end
            if ssdsop{i}(3)
                mdisloc.m(10,:)=1; % Set openning to 1
                for i2 = 1:size(mdisloc.m,2)
                    Uunit{i}{3}{i2} = disloc(mdisloc.m(:,i2),obs(1:2,:),nu); % displacements for unit slip
                end
            end
    end
end
                
%% Start core of inversion

sensitivityTest = 0; % Switch off sensitivity test at first iteration

while iKeep < invpar.nRuns  % While number of iterations is < than number of runs...
    if iKeep/invpar.TRuns == round(iKeep/invpar.TRuns) & iTemp < nTemp % Follow temperature schedule
        iTemp = iTemp + 1;
        T = invpar.TSchedule(iTemp);    % Assign temperature from T schedule
        if iKeep > 0
            model.trial = model.optimal;
        end
        if T ==1                        % Set Hyperparameter when T reaches 1 (Hyperparameter is currently not is use!!!)
            setHyperParameter = 1;
        else
            setHyperParameter = 0;
        end
    end
    
    
    if sum(iKeep == invpar.sensitivitySchedule)>0   % Check if it's time for sensitivity test based on schedule
        sensitivityTest = 1; % Switch on sensitivity test
    end
    
    %% Calculate 3D displacements from model
    
    UTot = zeros(3,nObs);   % Initialise matrix of modeled displacements (3 x number of observation points)
    Uinsar = cell(1,max(length(insar),length(gps)));
    for i = 1:invpar.nModels % For each source model...
        index1 = model.mIx(i);
        switch invpar.model{i}
            case 'MOGI'
                mFunc{i} = model.trial(index1:index1+3);    % Select source model parameters from all
                U = mogi(mFunc{i},obs,nu);               % Calculate 3D displacements
            case 'MCTG'
                mFunc{i} = model.trial(index1:index1+4);    % Select source model parameters from all
                U = mctigueSource(mFunc{i},obs(1:2,:),nu);        % Calculate 3D displacements
            case 'YANG'
                mFunc{i} = model.trial(index1:index1+7);    % Select source model parameters from all
                U = yangSource(mFunc{i},obs,nu);            % Calculate 3D displacements
            case 'PENN'
                mFunc{i}=model.trial(index1:index1+4);      % Select source model parameters from all
                U = pennySource(mFunc{i},obs,nu);           % Calculate 3D displacements
            case 'SILL'
                mFunc{i}=[model.trial(index1:index1+2); 0; model.trial(index1+3:index1+5); 0; 0; model.trial(index1+6)]; % Select source model parameters from all; Dip set to 0; Slip set to 0;
                U = disloc(mFunc{i},obs(1:2,:),nu); 
            case 'DIKE'
                mFunc{i}=[model.trial(index1:index1+6); 0; 0; model.trial(index1+7)]; % Select source model parameters from all; Slip set to 0;
                U = disloc(mFunc{i},obs(1:2,:),nu); 
            case 'FAUL'
                mFunc{i}=[model.trial(index1:index1+8);0]; % Select source model parameters from all; Opening set to 0;
                U = disloc(mFunc{i},obs(1:2,:),nu); 
            case 'DLOC'
                mFunc{i}=[model.trial(index1:index1+9)]; % Select source model parameters from all; Opening set to 0;
                U = disloc(mFunc{i},obs(1:2,:),nu); 
            case 'MDIK'
                U=zeros(3,size(obs,2));
                for i2=1:size(Uunit{i},2)
                    U = U + Uunit{i}{i2}*model.trial(index1+i2-1); 
                end
                mFunc{i}(10,:)=model.trial(index1:index1+i2-1);
             case 'MDLC'
                U=zeros(3,size(obs,2));
                Np=size(mFunc{i},2); % number of patches
                i_inc=0;
                for i2=1:3
                    if ssdsop{i}(i2)
                        for i3=1:Np
                            U = U + Uunit{i}{i2}{i3}*model.trial(index1+i_inc*Np+i3-1); 
                                                  
                        end
                        mFunc{i}(7+i2,:)=model.trial(index1+i_inc*Np:index1+(i_inc+1)*Np-1);
                        i_inc=i_inc+1;
                    end
                end
            case 'HING'
                mFunc{i}=model.trial(index1:index1+10);     % Select source model parameters from all
                U = hingedDikes(mFunc{i},obs(1:2,:),nu);    % Calculate 3D displacements
            case 'PCDM'
                mFunc{i}=model.trial(index1:index1+8); 
                % Calculate 3D displacements using Nikkhoo's code
                [ue,un,uv] = pCDM(obs(1,:),obs(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:), ...
                    mFunc{i}(4,:),mFunc{i}(5,:),mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),nu);
                U=[ue';un';uv'];
            case 'CDM'
                mFunc{i}=model.trial(index1:index1+9); 
                % Calculate 3D displacements using Nikkhoo's code
                [ue,un,uv] = CDM(obs(1,:),obs(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:), ...
                    mFunc{i}(4,:),mFunc{i}(5,:),mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),mFunc{i}(10,:),nu);
                U=[ue';un';uv'];
            case 'FECM'
                mu = modelInput.mu; % shear modulus
                lambda=2*mu*nu/(1-2*nu); % Lamé first constant
                mFunc{i}=model.trial(index1:index1+9);
                % Calculate 3D displacements using Nikkhoo's code
                [ue,un,uv,dV,DV,Ns] = fECM(xy(1,:),xy(2,:),mFunc{i}(1,:),mFunc{i}(2,:),mFunc{i}(3,:),mFunc{i}(4,:),mFunc{i}(5,:), ...
                    mFunc{i}(6,:),mFunc{i}(7,:),mFunc{i}(8,:),mFunc{i}(9,:),mFunc{i}(10,:),mu,lambda,'C');
                U=[ue';un';uv'];       
        end
        if isempty(invpar.modelinsarID{i})
            UTot = UTot + U; % Calculate total displacement from sum of displacement from each source
        else
            for i2=1:length(invpar.modelinsarID{i}) % only add for selected interferograms
                if isempty(Uinsar{invpar.modelinsarID{i}(i2)})
                    Uinsar{invpar.modelinsarID{i}(i2)}= U;
                else
                    Uinsar{invpar.modelinsarID{i}(i2)}=Uinsar{invpar.modelinsarID{i}(i2)} + U;
                end
            end
        end
    end
    
    insarParIx = model.mIx(end); % Identify first model parameter not related to source model (e.g., InSAR offset, ramp, etc.)
    
    
    %% Convert 3D displacement to LOS displacement and calculate residuals
    resExp = 0; % Initialise (Gm - d) * InvCov * (Gm - d)'
    
    if ~isempty(insar)
        
        ULos = []; % Initialise line-of-sight Gm vector
        
        for j = 1 : length(insar)
            U=UTot;
            if ~isempty(Uinsar{j})
                U=U+Uinsar{j};
            end
            UEast = -cosd(insar{j}.dHeading).* sind(insar{j}.dIncidence); % East unit vector
            UNorth = sind(insar{j}.dHeading).* sind(insar{j}.dIncidence); % North unit vector
            UVert = cosd(insar{j}.dIncidence); % Vertical unit vector
            
            ULos{j} = UEast.* U(1,insar{j}.ix) + ...
                UNorth.* U(2,insar{j}.ix) + ...             % Convert to line of sight displacement
                UVert.* U(3,insar{j}.ix);
            
            if insar{j}.constOffset == 'y'
                ULos{j} = ULos{j} + model.trial(insarParIx);  % Add constant offset
                
                insarParIx = insarParIx + 1; % Change model parameter index for next step
            end
            
            if insar{j}.rampFlag == 'y'
                ULos{j} = ULos{j} + model.trial(insarParIx)*insar{j}.obs(:,1)' + ...
                    model.trial(insarParIx+1)*insar{j}.obs(:,2)'; % Add linear ramp
                
                insarParIx = insarParIx + 2; % Change model parameter index for next step if necessary
            end
            
            resInsar{j} = (ULos{j} - insar{j}.dLos); % Calculate (Gm - d), residuals
            resExp = resExp + resInsar{j}* insar{j}.invCov* resInsar{j}'; % (Gm - d) * InvCov * (Gm - d)'
        end
    end
    
    
    %% Calculate GPS residuals
    
    if ~isempty(gps)
        for j = 1 : length(gps)
            if ~isempty(gps{j})
                U=UTot(1:3,gps{1}.ix);      
                if ~isempty(Uinsar{j})
                    U=U+Uinsar{j}(1:3,gps{1}.ix);
                end
                rGps = U - gps{j}.displacements(1:3,:); % Residual GPS displacement
                nonanix=~isnan(rGps(:));
                resExp = resExp + rGps(nonanix)' * gps{j}.invCov(nonanix,nonanix) * rGps(nonanix) * gps{j}.weight; % (Gm - d) * InvCov * (Gm - d)'
            end
        end
    end
    
    %% Calculate ENU residuals (treated like GPS-style point data)
    if exist('enu','var') && ~isempty(enu)
        for j = 1:length(enu)
            if ~isempty(enu{j})
                U = UTot(1:3, enu{1}.ix);
    
                % If you are using modelinsarID coupling, only add if appropriate.
                % Usually for GPS/ENU you should NOT index Uinsar by j (j is epoch),
                % so keep this disabled unless you know you want it.
                % if ~isempty(Uinsar{j})
                %     U = U + Uinsar{j}(1:3, enu{1}.ix);
                % end
    
                rEnu = U - enu{j}.displacements(1:3,:);     % 3 x nSites residuals
                nonanix = ~isnan(rEnu(:));
    
                % Weight (optional, match GPS behavior)
                if ~isfield(enu{j}, 'weight') || isempty(enu{j}.weight)
                    enu{j}.weight = 1;
                end
    
                resExp = resExp + rEnu(nonanix)' * ...
                    enu{j}.invCov(nonanix, nonanix) * ...
                    rEnu(nonanix) * enu{j}.weight;
            end
        end
    end

    %% Continue inversion ...
    
    if setHyperParameter == 1
        %hyperPrev = resExp/nObs; % set hyperparameter on first reaching T=1;
        hyperPrev = 1; % set hyperparameter to 1;
        model.trial(end) = log10(hyperPrev);
        setHyperParameter = 0; % Switch setHyperParameter off
    end
    
    if isempty(insar)
        hyperParam = 1;
    else
        hyperParam = 1;
        %hyperParam = 10^model.trial(end);
    end
    
    % !! Currently hyperparameter is set to 1
    P = -resExp/(2*hyperParam); % Probability is exp of P
    
    if iKeep>0
        PRatio = (hyperPrev/hyperParam)^(nObs/2)*exp((P-PPrev)/T);  % Probability ratio
    else
        PRatio=1; % Set to 1 for first iteration (always keep first iteration)
    end
    
    %% Perform sensitivity test if necessary and change step size
    
    if sensitivityTest > 1
        probSens(sensitivityTest-1) = PRatio; % Assign probability to current model parameter
        if sensitivityTest > nModel % Check if sensitivity test has finished
            if iKeepSave > 0
                rejectionRatio = (iReject - iRejectSave)/(iKeep - iKeepSave); % Calculate rejection rate
                probTarget = probTarget * rejectionRatio * 1/0.77; % Adjust target probability to reach 77% rejection rate
                probTarget(probTarget<1e-06) = 1e-06; % Prevent from reaching zero.
            end
            sensitivityTest = 0;    % Swtich off sensitivity test
            probSens(probSens > 1) = 1./probSens(probSens > 1);
            PDiff = probTarget - probSens;
            indexP = PDiff > 0; % Select model parameters for which to adjust model step
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/probTarget*2);  % assign new model step
            indexP = PDiff < 0; % Select remaining model parameters
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/(1-probTarget)*2); % assign new model step
            model.step(model.step > model.range) = model.range(model.step > model.range); % Check if step is within range
            iKeepSave = iKeep;
            iRejectSave = iReject;
        end
        
    else
        iKeep = iKeep + 1;
        if PRatio >= rand(1,1)  % If condions are met, keep model trial
            model.m = model.trial; % Substitute m with model trial
            mKeep(:,iKeep) = model.m;   % Keep model trial
            PKeep(:,iKeep) = P;         % P of current model
            
            PPrev = P;  % Assign current P to PPrev for next trial
            hyperPrev = hyperParam; % Assign current Hyperparameter for next trial
            
            if P > POpt   % Update current optimal model if likelihood is higher
                model.optimal = model.m;
                model.funcOpt = mFunc;
                POpt = P;
            end
        else                    % Reject model trial and keep previous model
            iReject = iReject + 1;
            mKeep(:,iKeep) = mKeep(:,iKeep-1);
            PKeep(:,iKeep) = PKeep(:,iKeep-1);
        end
        
        if iKeep/invpar.nSave == round(iKeep/invpar.nSave) % display and save results at regular intervals (1000 or 10000 iterations)
            if iKeep >= 20000           % Increase time step for saving/displaying after 20000 iterations
                invpar.nSave = 10000;
            end
            
            % Print current status of inversion to screen
            disp('=========================================================')
            disp(['Model: ',invpar.model{:}])
            disp([num2str(iKeep),' model trials. Optimal Prob = exp(',num2str(POpt),')'])
            disp(['Hyperparameter=',num2str(hyperParam)])
            disp([num2str(iReject),' models rejected:', num2str((iReject/iKeep)*100),'% of model trials.'])
            
            % allocate space for next blocks to keep
            mKeep(:,iKeep + invpar.nSave) = 0;
            PKeep(:,iKeep + invpar.nSave) = 0;
            
            % Save results to temporary file for insepction during
            % inversion
            save([outputDir,'/temporary.mat'], 'geo', 'mKeep', 'PKeep', 'model', 'gps', 'enu', 'insar', 'invpar', 'geo', 'modelInput','-v7.3');
            
            % Display current optimal model parameters on screen
            for i=1:length(invpar.model)
                if invpar.model{i} == 'MOGI'
                    fprintf('\nMOGI center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('MOGI center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('MOGI depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('MOGI volume change: %f\n',(model.funcOpt{i}(4,:)));
                elseif invpar.model{i} == 'YANG'
                    fprintf('\nYANG centroid X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('YANG centroid Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('YANG centroid depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('YANG major axis: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('YANG minor axis: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('YANG majax strike: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('YANG majax plunge: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('YANG DP/mu: %f\n',(model.funcOpt{i}(8,:)));
                elseif invpar.model{i} == 'MCTG'
                    fprintf('\nMCTG center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('MCTG center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('MCTG center depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('MCTG radius: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('MCTG DP/mu: %f\n',(model.funcOpt{i}(5,:)));
                elseif invpar.model{i} == 'PENN'
                    fprintf('\nPENNY center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('PENNY center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('PENNY center depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('PENNY radius: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('PENNY DP/mu: %f\n',(model.funcOpt{i}(5,:)));
                elseif invpar.model{i} == 'SILL'
                    fprintf('\nSILL length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('SILL width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('SILL depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('SILL strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('SILL edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('SILL edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('SILL opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif invpar.model{i} == 'DIKE'
                    fprintf('\nDIKE length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('DIKE width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('DIKE depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('DIKE dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('DIKE strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('DIKE edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('DIKE edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('DIKE opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif invpar.model{i} == 'FAUL'
                    fprintf('\nFAULT length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('FAULT width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('FAULT depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('FAULT dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('FAULT strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('FAULT edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FAULT edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FAULT strike-slip component: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('FAULT dip-slip component: %f\n',(model.funcOpt{i}(9,:)));
                elseif invpar.model{i} == 'DLOC'
                    fprintf('\nDISLOCATION length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('DISLOCATION width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('DISLOCATION depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('DISLOCATION dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('DISLOCATION strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('DISLOCATION edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('DISLOCATION edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('DISLOCATION strike-slip component: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('DISLOCATION dip-slip component: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('DISLOCATION opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif invpar.model{i} == 'MDIK'
                    fprintf('\nMULTI-PATCH DIKE Opening: %f to %f\n',min(model.funcOpt{i}(10,:)),max(model.funcOpt{i}(10,:)));    
                elseif invpar.model{i} == 'MDLC'
                    fprintf('\nMULTI-PATCH Dislocation Strike-slip: %f to %f\n',min(model.funcOpt{i}(8,:)),max(model.funcOpt{i}(8,:)));  
                    fprintf('MULTI-PATCH Dislocation Dip-slip: %f to %f\n',min(model.funcOpt{i}(9,:)),max(model.funcOpt{i}(9,:)));  
                    fprintf('MULTI-PATCH Dislocation Opening: %f to %f\n',min(model.funcOpt{i}(10,:)),max(model.funcOpt{i}(10,:)));    
                elseif invpar.model{i} == 'HING'
                    fprintf('\nDIKE1 edge center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('DIKE1 edge center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('DIKE1 length: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('DIKE1 width: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('DIKE1 depth: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('DIKE1 dip: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('\nDIKE1 opening: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('DIKE2 width: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('DIKE2 dip: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('DIKE2 opening: %f\n',(model.funcOpt{i}(10,:)));
                    fprintf('Strike: %f\n\n',(model.funcOpt{i}(11,:)));
                elseif strcmp(invpar.model{i},'PCDM')
                    fprintf('\nPCDM centroid X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('PCDM centroid Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('PCDM centroid depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('PCDM x rotation angle: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('PCDM y rotation angle: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('PCDM z rotation angle: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('PCDM x Potencies of the PTDs: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('PCDM y Potencies of the PTDs: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('PCDM z Potencies of the PTDs: %f\n',(model.funcOpt{i}(9,:)));
                elseif strcmp(invpar.model{i},'CDM')
                    fprintf('\nCDM centroid X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('CDM centroid Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('CDM centroid depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('CDM x rotation angle: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('CDM y rotation angle: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('CDM z rotation angle: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('CDM semi-axis x: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('CDM semi-axis y: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('CDM semi-axis z: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('CDM opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif strcmp(invpar.model{i},'FECM')
                    fprintf('\nFECM centroid X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('FECM centroid Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('FECM centroid depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('FECM x rotation angle: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('FECM y rotation angle: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('FECM z rotation angle: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FECM semi-axis x: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FECM semi-axis y: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('FECM semi-axis z: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('FECM pressure: %f\n',(model.funcOpt{i}(10,:)));
                end
            end
        end
    end
    
    if sensitivityTest > 0  % Perform sensitivity test (no models are kept during this phase!)
        randomStep = zeros(nModel,1);
        randomStep(sensitivityTest) = model.step(sensitivityTest) * sign(randn(1,1))/2; % Assign random step
        model.trial = model.m + randomStep; % New model trial
        % Check that new model trial is withing bounds
        if model.trial(sensitivityTest) > model.upper(sensitivityTest)
            model.trial(sensitivityTest) = model.trial(sensitivityTest) - model.step(sensitivityTest);
        end
        
        hyperParam = hyperPrev;
        sensitivityTest = sensitivityTest + 1; % Move index to that of next parameter until all parameters are done
    else
        randomStep = model.step.*(rand(nModel,1)-0.5)*2;     % Make random step
        model.trial = model.m + randomStep;                 % Assign new model trial to previous + random step
        % Check that new model trial is withing bounds
        model.trial(model.trial > model.upper) = 2 * model.upper(model.trial > model.upper) - ...
            model.trial(model.trial > model.upper);
        
        model.trial(model.trial < model.lower) = 2 * model.lower(model.trial < model.lower) - ...
            model.trial(model.trial < model.lower);
    end
end


%% Clean up and prepare results
mKeep(:, end - invpar.nSave) = []; % Remove unused preallocated memory
PKeep(:, end - invpar.nSave) = []; % Remove unused preallocated memory

results.mKeep = mKeep;
results.PKeep = PKeep;
results.model = model;
results.optimalmodel = model.funcOpt;

