function model = prepareModel(modelInput, invpar, insar, gps, saveName, restartFlag)

% Function to prepare model parameters
%
% Usage: model = prepareModel(modelInput, invpar, insar, gps)
% Input Parameters:
%       modelInput: parameters read from input file
%       invpar: inversion parameters
%
% Output Parameters:
%       model: structure containing model settings to be used for inversion
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
% Last update: 16 Jan, 2023

global outputDir  % Set global variables

%% Initialize model vectors
mIx = zeros(invpar.nModels+1, 1);
mIx(1) = 1;
model.m = zeros(500,1);
model.step = model.m;
model.lower = model.m;
model.upper = model.m;

%% Assign model parameters from input file


i_mogi=0;
i_yang=0;
i_mctigue=0;
i_sill=0;
i_dike=0;
i_penny=0;
i_fault=0;
i_disloc=0;
i_mdike=0;
i_mdisloc=0;
i_hing=0;
i_pcdm=0;
i_cdm=0;
i_fecm=0;
  for i = 1 : invpar.nModels
    index1 = mIx(i);
    switch invpar.model{i}
        case 'MOGI'
            i_mogi=i_mogi+1;
            nParameters = 4;
            if isstruct(modelInput.mogi) % in case no numbers assigned in input file (old format)
                atemp=modelInput.mogi;
                modelInput.mogi=[];
                modelInput.mogi{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.mogi{i_mogi}.start;
            model.step(index1:index2) = modelInput.mogi{i_mogi}.step;
            model.lower(index1:index2) = modelInput.mogi{i_mogi}.lower;
            model.upper(index1:index2) = modelInput.mogi{i_mogi}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['MOGI',num2str(i_mogi)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'DV'};
        case 'YANG'
            i_yang=i_yang+1;
            nParameters = 8;
            if isstruct(modelInput.yang) % in case no numbers assigned in input file (old format)
                atemp=modelInput.yang;
                modelInput.yang=[];
                modelInput.yang{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.yang{i_yang}.start;
            model.step(index1:index2) = modelInput.yang{i_yang}.step;
            model.lower(index1:index2) = modelInput.yang{i_yang}.lower;
            model.upper(index1:index2) = modelInput.yang{i_yang}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Yang',num2str(i_yang)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'majAx'; 'a/r'; ...
                'strike'; 'Plunge'; 'DP/mu'};
        case 'MCTG'
            i_mctigue=i_mctigue+1;
            nParameters = 5;
            if isstruct(modelInput.mctigue) % in case no numbers assigned in input file (old format)
                atemp=modelInput.mctigue;
                modelInput.mctigue=[];
                modelInput.mctigue{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.mctigue{i_mctigue}.start;
            model.step(index1:index2) = modelInput.mctigue{i_mctigue}.step;
            model.lower(index1:index2) = modelInput.mctigue{i_mctigue}.lower;
            model.upper(index1:index2) = modelInput.mctigue{i_mctigue}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Mctigue',num2str(i_mctigue)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'Radius'; 'DP/mu'};
        case 'PENN'
            i_penny=i_penny+1;
            nParameters = 5;
            if isstruct(modelInput.penny) % in case no numbers assigned in input file (old format)
                atemp=modelInput.penny;
                modelInput.penny=[];
                modelInput.penny{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.penny{i_penny}.start;
            model.step(index1:index2) = modelInput.penny{i_penny}.step;
            model.lower(index1:index2) = modelInput.penny{i_penny}.lower;
            model.upper(index1:index2) = modelInput.penny{i_penny}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1); 
            model.modelName(index1:index2)={['Penny',num2str(i_penny)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'Radius'; 'DP/mu'};
        case 'SILL'
            i_sill=i_sill+1;
            nParameters = 7;
            if isstruct(modelInput.sill) % in case no numbers assigned in input file (old format)
                atemp=modelInput.sill;
                modelInput.sill=[];
                modelInput.sill{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.sill{i_sill}.start;
            model.step(index1:index2) = modelInput.sill{i_sill}.step;
            model.lower(index1:index2) = modelInput.sill{i_sill}.lower;
            model.upper(index1:index2) = modelInput.sill{i_sill}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Sill',num2str(i_fault)]};
            model.parName(index1:index2) = {'Length'; 'Width'; 'Depth'; 'Strike'; 'X'; 'Y'; 'Opening'};
        case 'DIKE'
            i_dike=i_dike+1;
            nParameters = 8;
            if isstruct(modelInput.dike) % in case no numbers assigned in input file (old format)
                atemp=modelInput.dike;
                modelInput.dike=[];
                modelInput.dike{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.dike{i_dike}.start;
            model.step(index1:index2) = modelInput.dike{i_dike}.step;
            model.lower(index1:index2) = modelInput.dike{i_dike}.lower;
            model.upper(index1:index2) = modelInput.dike{i_dike}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Dike',num2str(i_dike)]};
            model.parName(index1:index2) = {'Length'; 'Width'; 'Depth'; 'Dip'; 'Strike'; 'X'; 'Y'; 'Opening'};
        case 'FAUL'
            i_fault=i_fault+1;
            nParameters = 9;
            if isstruct(modelInput.fault) % in case no numbers assigned in input file (old format)
                atemp=modelInput.fault;
                modelInput.fault=[];
                modelInput.fault{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.fault{i_fault}.start;
            model.step(index1:index2) = modelInput.fault{i_fault}.step;
            model.lower(index1:index2) = modelInput.fault{i_fault}.lower;
            model.upper(index1:index2) = modelInput.fault{i_fault}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Fault',num2str(i_fault)]};
            model.parName(index1:index2) = {'Length'; 'Width'; 'Depth'; 'Dip'; 'Strike'; 'X'; 'Y'; 'StrSlip'; 'DipSlip'};
        case 'DLOC'
            i_disloc=i_disloc+1;
            nParameters = 10;
            if isstruct(modelInput.disloc) % in case no numbers assigned in input file (old format)
                atemp=modelInput.disloc;
                modelInput.disloc=[];
                modelInput.disloc{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.disloc{i_disloc}.start;
            model.step(index1:index2) = modelInput.disloc{i_disloc}.step;
            model.lower(index1:index2) = modelInput.disloc{i_disloc}.lower;
            model.upper(index1:index2) = modelInput.disloc{i_disloc}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Disloc',num2str(i_disloc)]};
            model.parName(index1:index2) = {'Length'; 'Width'; 'Depth'; 'Dip'; 'Strike'; 'X'; 'Y'; 'StrSlip'; 'DipSlip'; 'Opening'};
        case 'MDIK'
            i_mdike=i_mdike+1;
            if isstruct(modelInput.mdike) % in case no numbers assigned in input file (old format)
                atemp=modelInput.mdike;
                modelInput.mdike=[];
                modelInput.mdike{1}=atemp;
            end
            mdike=load(modelInput.mdike{i_mdike}.patchfile);
            nParameters=size(mdike.m,2);
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.mdike{i_mdike}.start;
            model.step(index1:index2) = modelInput.mdike{i_mdike}.step;
            model.lower(index1:index2) = modelInput.mdike{i_mdike}.lower;
            model.upper(index1:index2) = modelInput.mdike{i_mdike}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['MDike',num2str(i_mdike)]};
            model.parName(index1:index2) = {'Opening'}; 
        case 'MDLC'
            i_mdisloc=i_mdisloc+1;
            if isstruct(modelInput.mdisloc) % in case no numbers assigned in input file (old format)
                atemp=modelInput.mdisloc;
                modelInput.mdisloc=[];
                modelInput.mdisloc{1}=atemp;
            end
            mdisloc=load(modelInput.mdisloc{i_mdisloc}.patchfile);
            ssdsop=modelInput.mdisloc{i_mdisloc}.ss_ds_op;
            nParameters=size(mdisloc.m,2)*sum(ssdsop);
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.mdisloc{i_mdisloc}.start;
            model.step(index1:index2) = modelInput.mdisloc{i_mdisloc}.step;
            model.lower(index1:index2) = modelInput.mdisloc{i_mdisloc}.lower;
            model.upper(index1:index2) = modelInput.mdisloc{i_mdisloc}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['MDisloc',num2str(i_mdisloc)]};
            model.parName(index1:index2) = {''}; 
        case 'HING'
            i_hing=i_hing+1;
            nParameters = 11;
            if isstruct(modelInput.hing) % in case no numbers assigned in input file (old format)
                atemp=modelInput.hing;
                modelInput.hing=[];
                modelInput.hing{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.hing{i_hing}.start;
            model.step(index1:index2) = modelInput.hing{i_hing}.step;
            model.lower(index1:index2) = modelInput.hing{i_hing}.lower;
            model.upper(index1:index2) = modelInput.hing{i_hing}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['Hinge',num2str(i_hing)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Length'; 'Width1'; 'Depth1'; 'Phi1'; 'Open1'; 'Width2'; 'Phi2';'Open2'; 'Strike'};
        case 'PCDM'
            i_pcdm=i_pcdm+1;
            nParameters = 9;
            if isstruct(modelInput.pcdm) % in case no numbers assigned in input file (old format)
                atemp=modelInput.pcdm;
                modelInput.pcdm=[];
                modelInput.pcdm{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.pcdm{i_pcdm}.start;
            model.step(index1:index2) = modelInput.pcdm{i_pcdm}.step;
            model.lower(index1:index2) = modelInput.pcdm{i_pcdm}.lower;
            model.upper(index1:index2) = modelInput.pcdm{i_pcdm}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['PCDM',num2str(i_pcdm)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'OmegaX'; 'OmegaY'; 'OmegaZ'; 'DV_X'; 'DV_Y'; 'DV_Z'};
        case 'CDM'
            i_cdm=i_cdm+1;
            nParameters = 10;
            if isstruct(modelInput.cdm) % in case no numbers assigned in input file (old format)
                atemp=modelInput.cdm;
                modelInput.cdm=[];
                modelInput.cdm{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.cdm{i_cdm}.start;
            model.step(index1:index2) = modelInput.cdm{i_cdm}.step;
            model.lower(index1:index2) = modelInput.cdm{i_cdm}.lower;
            model.upper(index1:index2) = modelInput.cdm{i_cdm}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['CDM',num2str(i_cdm)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'Rotation X'; ...
                'Rotation Y'; 'Rotation Z'; 'Semi-axis X'; 'Semi-axis Y'; ...
                'Semi-axis Z'; 'Opening'};
        case 'FECM'
            i_fecm=i_fecm+1;
            nParameters = 10;
            if isstruct(modelInput.fecm) % in case no numbers assigned in input file (old format)
                atemp=modelInput.fecm;
                modelInput.fecm=[];
                modelInput.fecm{1}=atemp;
            end
            index2 = index1 + nParameters - 1;
            model.m(index1:index2) = modelInput.fecm{i_fecm}.start;
            model.step(index1:index2) = modelInput.fecm{i_fecm}.step;
            model.lower(index1:index2) = modelInput.fecm{i_fecm}.lower;
            model.upper(index1:index2) = modelInput.fecm{i_fecm}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.modelName(index1:index2)={['FECM',num2str(i_fecm)]};
            model.parName(index1:index2) = {'X'; 'Y'; 'Depth'; 'Rotation X'; ...
                'Rotation Y'; 'Rotation Z'; 'Semi-axis X'; 'Semi-axis Y'; ...
                'Semi-axis Z'; 'Pressure'};
        otherwise
            error('Invalid model')
    end
    mIx(i+1) = mIx(i) + nParameters;
  end

if strcmpi(restartFlag,'y') % use saved values from previous run
    filePath = fullfile(pwd, outputDir, saveName);
    oldModel = load(filePath, 'invResults');
    display(filePath);
    model.m(1:index2)=oldModel.invResults.model.m(1:index2);
end

% Add other parameters to invert for (e.g., InSAR constant offset, ramp,
% etc.)
clear index1
index1 = index2+1;
clear index2

nParameters = 0;
insarParName = {};

for i=1:length(insar)
    
    if insar{i}.constOffset == 'y';
        nParameters = nParameters + 1; % Constant offset
        if i == 1
            insarParName = {'Constant'};
        else
            insarParName = [insarParName, 'Constant'];
        end
    end
    
    if insar{i}.rampFlag == 'y'
        nParameters = nParameters + 2;  % Linear ramp
        insarParName = [insarParName, 'X-ramp', 'Y-ramp';];
    end
end

index2 = index1 + nParameters; % this also adds hyperparameter at end of each model vector

model.m(index1:index2) = zeros(nParameters+1,1);
model.step(index1:index2) = ones(nParameters+1,1)*1e-3;
model.lower(index1:index2) = ones(nParameters+1,1)*-5e-1;
model.upper(index1:index2) = ones(nParameters+1,1)*5e-1;
if ~isempty(insar)
    model.modelName(index1:index2-1)={'InSAR'};
    model.parName(index1:index2-1) = insarParName;
end

% Clear unused values
model.m = model.m(1:index2);
model.step = model.step(1:index2);
model.lower = model.lower(1:index2);
model.upper = model.upper(1:index2);

model.mIx = mIx;
