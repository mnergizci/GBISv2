function results = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs)
%%calculation of the resolution of patches and illustrate in GBIS
%MNergizci August 2024. s
global outputDir  % Set global variables
figdir = [outputDir, '/Figures/'];
nu = modelInput.nu; % Poisson's ratio
i_mdisloc=0;
for i = 1:invpar.nModels % For each source model...
    index1 = model.mIx(i);
    switch invpar.model{i}
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

%% Generate G and resolution matrix
% Initialize los_vector_inSAR as an empty matrix
los_vector_inSAR = [];
invCov_total = [];
% Generate los_vector_inSAR
for j = 1:length(insar)
    UEast = cosd(insar{j}.dHeading) .* sind(insar{j}.dIncidence); % East unit vector
    UNorth = sind(insar{j}.dHeading) .* sind(insar{j}.dIncidence); % North unit vector
    UVert = cosd(insar{j}.dIncidence); % Vertical unit vector
    % Combine the vectors into a 3xj matrix and store in los_vector_inSAR
    los_vector_inSAR = [los_vector_inSAR, [UEast; UNorth; UVert]];
    invCov_total=blkdiag(invCov_total, insar{j}.invCov);
end

% Initialize G_total as a cell array
G_total = cell(1, 3);

% Loop through models and calculate G matrix
for i = 1:invpar.nModels
    switch invpar.model{i}
        case 'MDLC'
            disp(' ')
            disp('resolution matrix is calculated for MDLC!..')
            Np = size(mFunc{i}, 2); % number of patches
            
            for i2 = 1:3 % iteration for strike and dip slip
                if ssdsop{i}(i2)
                    G_total{i2} = cell(1, Np); % Initialize cells for G_total
                    for i3 = 1:Np % iteration over number of patches
                        U_temp = Uunit{i}{i2}{i3}; % 3*nObs, 3 is ENU component
                        % Compute the dot product of U_temp and los_vector_inSAR
                        G_total{i2}{i3} = sum(U_temp .* los_vector_inSAR, 1); % Summing along the columns
                    end
                end
            end
    end
end
% Convert G_total{1} to a matrix G_ss
if ~isempty(G_total{1})
    G_ss = cell2mat(G_total{1}');
    G_ss_t= G_ss';
else
    G_ss = [];
end
% Convert G_total{2} to a matrix G_ds
if ~isempty(G_total{2})
    G_ds = cell2mat(G_total{2}');
    G_ds_t= G_ds';
else
    G_ds = [];
end

% %rake value  not meaningful it gives you G_ss_t essentially
% rake = zeros(1,Np);
% G=G_ss_t*diag(cosd(rake))+G_ds_t*diag(sind(rake));

G=G_ss_t;
GT=G';
Gg=inv(GT*invCov_total*G)*GT*invCov_total;
resolution_matrices=Gg*G;
patch_resolution = diag(resolution_matrices);

%% plotting
climres= [0,1];
eqs_sw = 'n';
chot = hot;
chot(1:40,:) = [];
chot = [chot;ones(30,3)];

view_angle = [-4, 25];
shiftorigin = [0, 0];

m_res=mdisloc.m;
m_res(10,:)=patch_resolution;
setplotattr(m_res, chot, climres, eqs_sw, view_angle, shiftorigin);
print('-dpng', [figdir, 'resolution_matrix', num2str(i), '.png'], '-r900');
exportgraphics(gcf, [figdir, 'resolution_matrix', num2str(i), '.pdf'], 'ContentType', 'vector');

%saving the results
results=patch_resolution;



function setplotattr(m_plot, chot, clim, eqs_sw, view_angle, shiftorigin)
    colormap(chot);
    colormap(flipud(colormap));
    H = drawmodel(m_plot, 'color', [0.8, 0.8, 0.8], 'updipline', 'no', 'openingcolor', clim, 'outline', 'yes');
    set(H, 'linewidth', 0.15, 'facealpha', 1);

    if strcmpi(eqs_sw, 'y')
        [xy_patch, depth] = plot_eqs('n', 0, 0, 'y');
        plot3(xy_patch(:, 1) / 1000 + shiftorigin(1), xy_patch(:, 2) / 1000 + shiftorigin(2), -depth / 1000, 'k.', 'markersize', 5);
    end
    view(view_angle);

    if view_angle(1) == 0
        set(gca, 'DataAspectRatio', [1, 1, 1.7]);
        set(gca, 'xtick', [0:sind(m_plot(5)):8.5]);
        set(gca, 'xticklabel', strsplit(num2str([0:14])));
        xlabel('Along strike (km)');
    else
        axis equal;
        xlabel('Easting (km)');
    end

    ylabel('Northing (km)');
    zlabel('Depth (km)');
    grid on;
    box on;
    ztick = get(gca, 'zticklabel');
    for i1 = 1:length(ztick)
        ztick{i1} = num2str(abs(str2num(ztick{i1})));
    end
    set(gca, 'zticklabel', ztick);

    c = colorbar('vertical');
    if view_angle(1) == 0
        set(c, 'position', [0.930, 0.3800, 0.02900, 0.3000]);
    else
        set(c, 'position', [0.903 0.398 0.013 0.422]);
    end
%     xlabel(c, 'Opening (m)');
