function plot_grind(inputfile, dirname, burn_in, eqs_sw, flat_sw)
    % PLOT_GRIND

    if nargin < 1
        % inputfile='invert_1_2_3_GPS_Z_M';
        inputfile = 'invert_3_6_9_Z';
    end

    if nargin < 2
        % dirname=['dike_fixed'];
        dirname = ['/nobackup/mnergizci/GBIS_V2.0.1/Example/sparse_TR_exercise/'];
    end

    if nargin < 3
        burn_in = 200000;
    end
    
    if nargin < 4
        eqs_sw = 'n';
    end

    if nargin < 5
        flat_sw = 'n';
    end

    ldirname = [dirname, '/', inputfile, '/'];
    load([ldirname, inputfile, '.mat']);

    N = round(sum(invResults.PKeep ~= 0) / 2);

      % Number of empty cells at the end of mKeep and pKeep
    if invpar.nRuns < 10000
        blankCells = 999; 
    else
        blankCells = 9999;
    end  
    
    
    chot = hot(64);
    chot = flipud(chot(1:end, :));
    chot = chot(1:end-10, :);
    climOP = [0, 10];
    climSS = [0, 10];
    climDS = [-2, 2];
    climstd = [0, 5];
    shiftorigin = [0, 0];

    figdir = [ldirname, 'Figures/'];
    close all;

    for i = 1:length(invpar.model)
        if strcmpi(invpar.model{i}, 'MDLC')
            m_plot = invResults.optimalmodel{i};
            m_plot([1:3, 6, 7], :) = m_plot([1:3, 6, 7], :) / 1000; % changing to meter to km for illustration!
            m_plot(6, :) = m_plot(6, :) + shiftorigin(1); % reference to end of dike
            m_plot(7, :) = m_plot(7, :) + shiftorigin(2);
            ss = sum(m_plot(8, :)) ~= 0; % checking the which type of the slip is modelled and unmodelled during inversion.
            ds = sum(m_plot(9, :)) ~= 0; % return 1 or 0 as logic.
            op = sum(m_plot(10, :)) ~= 0;
            n_patch = size(m_plot, 2);
            m_keep_orig = invResults.mKeep(:, invResults.PKeep ~= 0);
%             m_keep_x = (m_keep(model.mIx(i):model.mIx(i+1)-1, end-N+1:end));
            m_keep = (m_keep_orig(model.mIx(i):model.mIx(i+1)-1, burn_in:end-blankCells));

            if abs(m_plot(4, 1)) < 10
                view_angle = [-4, 25];
            else
                view_angle = [0, 120];
                view_angle = [-4, 25];
            end
            if strcmpi(flat_sw, 'y')
                view_angle = [0, 0];
            end

            if ss
                %% MAP
                figure
                m_plot(10, :) = abs(m_plot(8, :));
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'OptimalSS', num2str(i), '.png'], '-r900');
%                 exportgraphics(gcf, [figdir, 'OptimalSS', num2str(i), '.pdf'], 'ContentType', 'vector');

                %% Median
                ix = 1;
                m_plot(10, :) = abs(median(m_keep(ix:ix+n_patch-1, :), 2));
                figure
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'MedianSS', num2str(i), '.png'], '-r900');
%                 exportgraphics(gcf, [figdir, 'MedianSS', num2str(i), '.pdf'], 'ContentType', 'vector');


%                 %% Bootstrap for Confidence Interval of the Median
%                 n_resamples = 1000;  % Number of bootstrap resamples
%                 alpha = 0.05;  % Significance level for the 95% confidence interval
%                 boot_medians = zeros(n_resamples, n_patch);
%                 
%                 % Perform bootstrapping
%                 for b = 1:n_resamples
%                     % Resample with replacement
%                     resample_ix = randi(size(m_keep(ix:ix+n_patch-1, :), 2), 1, size(m_keep(ix:ix+n_patch-1, :), 2));
%                     resample_data = m_keep(ix:ix+n_patch-1, resample_ix);
%                     
%                     % Calculate the median of the resampled data
%                     boot_medians(b, :) = abs(median(resample_data, 2));
%                 end
%                 
%                 % Calculate confidence intervals for each patch
%                 lower_bound = prctile(boot_medians, 100 * (alpha / 2), 1);  % Lower bound (2.5 percentile for 95% CI)
%                 upper_bound = prctile(boot_medians, 100 * (1 - alpha / 2), 1);  % Upper bound (97.5 percentile for 95% CI)
%                 boots_median = abs(median(m_keep(ix:ix+n_patch-1, :), 2));  % Original median
%                 
%                 % Plotting the bootstrap median
%                 figure
%                 m_plot(10, :) = boots_median;  % Use the original median for visualization
%                 setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
%                 print('-dpng', [figdir, 'bootsmedianss', num2str(i), '.png'], '-r900');

                %% Mean
                ix = 1;
                m_plot(10, :) = abs(mean(m_keep(ix:ix+n_patch-1, :), 2));
                figure
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'MeanSS', num2str(i), '.png'], '-r900');
%                 exportgraphics(gcf, [figdir, 'MeanSS', num2str(i), '.pdf'], 'ContentType', 'vector');

                %% std
                ix = 1;
                m_plot(10, :) = abs(std(m_keep(ix:ix+n_patch-1, :), 0, 2));
                figure
                setplotattr(m_plot, chot, climstd, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'stdSS', num2str(i), '.png'], '-r900');
%                 exportgraphics(gcf, [figdir, 'stdSS', num2str(i), '.pdf'], 'ContentType', 'vector');
            end

            if ds
                %% MAP
                figure
                m_plot(10, :) = (m_plot(9, :));
                setplotattr(m_plot, chot, climDS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'OptimalDS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'OptimalDS', num2str(i), '.pdf'], 'ContentType', 'vector');

                %% Median
                ix = n_patch + 1;
                m_plot(10, :) = (median(m_keep(ix:ix+n_patch-1, :), 2));
                figure
                setplotattr(m_plot, chot, climDS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'MedianDS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'MedianDS', num2str(i), '.pdf'], 'ContentType', 'vector');
                
                %% Mean
                ix = n_patch + 1;
                m_plot(10, :) = (mean(m_keep(ix:ix+n_patch-1, :), 2));
                figure
                setplotattr(m_plot, chot, climDS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'MeanDS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'MeanDS', num2str(i), '.pdf'], 'ContentType', 'vector');

                %% std
                ix = n_patch + 1;
                m_plot(10, :) = (std(m_keep(ix:ix+n_patch-1, :), 0, 2));
                figure
                setplotattr(m_plot, chot, climstd, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'stdDS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'stdDS', num2str(i), '.pdf'], 'ContentType', 'vector');
            end

            if ds && ss
                %% MAP
                figure
                m_plot(10, :) = sqrt((abs(m_plot(8, :)).^2) + (abs(m_plot(9, :)).^2));
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'Optimal_totalS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'Optimal_totalS', num2str(i), '.pdf'], 'ContentType', 'vector');

                %% Median
                ix_ss = 1; % Start index for ss
                ix_ds = n_patch + 1; % Start index for ds
                median_ss = abs(median(m_keep(ix_ss:ix_ss+n_patch-1, :), 2));
                median_ds = abs(median(m_keep(ix_ds:ix_ds+n_patch-1, :), 2));
                median_total = sqrt(median_ss.^2 + median_ds.^2);
                m_plot(10, :) = median_total;
                figure
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'Median_totalS', num2str(i), '.png']);
%                 exportgraphics(gcf, [figdir, 'Median_totalS', num2str(i), '.pdf'], 'ContentType', 'vector');
                
                %% Mean
                ix_ss = 1; % Start index for ss
                ix_ds = n_patch + 1; % Start index for ds
                mean_ss = abs(mean(m_keep(ix_ss:ix_ss+n_patch-1, :), 2));
                mean_ds = abs(mean(m_keep(ix_ds:ix_ds+n_patch-1, :), 2));
                m_plot(10, :) = sqrt(mean_ss.^2 + mean_ds.^2);
                figure
                setplotattr(m_plot, chot, climSS, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'Mean_totalS', num2str(i), '.png']);

                %% std
                ix_ss = 1; % Start index for ss
                ix_ds = n_patch + 1; % Start index for ds
                % Calculate standard deviations
                std_ss = abs(std(m_keep(ix_ss:ix_ss+n_patch-1, :), 0, 2));
                std_ds = abs(std(m_keep(ix_ds:ix_ds+n_patch-1, :), 0, 2));
                m_plot(10, :) = sqrt(std_ss.^2 + std_ds.^2);
                figure
                setplotattr(m_plot, chot, climstd, eqs_sw, view_angle, shiftorigin);
                print('-dpng', [figdir, 'std_totalS', num2str(i), '.png']);
% %                 exportgraphics(gcf, [figdir, 'std_totalS', num2str(i), '.pdf'], 'ContentType', 'vector');
                
                %% marginal pdf
                fontsize_plot = 16;
                set(0,'DefaultAxesFontSize',fontsize_plot)
                set(0,'defaultAxesFontName', 'Times New Roman')
                disp('calculation of random patch pdfs')
                %you can illustrate the 5 index which includes 5 values.
                [~, sortIndex]=sort(median_ss, 'descend');
                maxIndex=sortIndex(1:5);
                slip_for_marginal_PDF_plotting=abs(m_keep_orig(maxIndex,:));
                
                %select the 10 random index
                randomIndex = randperm(length(median_ss), 10);
                slip_for_marginal_PDF_plotting=abs(m_keep(randomIndex,:));    

                if size(slip_for_marginal_PDF_plotting, 2) > size(slip_for_marginal_PDF_plotting, 1)   % transpose if necessary
                    slip_for_marginal_PDF_plotting = slip_for_marginal_PDF_plotting';
                end
                hFig = figure('position', [100, 300, 1600, 1200]);
                [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower(slip_for_marginal_PDF_plotting,'plot_color'); %%thanks to slipBERI (Ruth Amey and David Bakeart).
                
                for i = 1: size(slip_for_marginal_PDF_plotting,2)
                   strings{i} = (['P', num2str(sortIndex(i))]) ;
                   strings_units{i} = (['P', num2str(sortIndex(i)), ' slip (m) ']);                      % (['Patch ', num2str(i), '  ';'slip [m] '])
                end
                
                for k = 1 : (size(slip_for_marginal_PDF_plotting,2)-1)                      % loop to (n-1) because we want to put the labels on the LHS from slip patch 2:end, not 1:end. Then add on the bottom label for last slip patch.
                   set(get(PAx(k),'xlabel'),'string',strings_units{k},'fontsize',fontsize_plot-6);         %set(get(PAx(k),'xlabel'),'string',strings_units{k},'fontsize',fontsize_plot-1)
                   set(get(AX(k,1),'ylabel'),'string',strings{k+1},'fontsize',fontsize_plot-6);     % set(get(AX(k,1),'ylabel'),'string',strings{k+1},'fontsize',fontsize_plot-1)
                end
                set(get(PAx( size(slip_for_marginal_PDF_plotting,2)),'xlabel'),'string',strings_units{size(slip_for_marginal_PDF_plotting,2)},'fontsize',fontsize_plot-6);      % add on bottom label for last slip patch, name is the last value in the matrix 'strings'
                print(hFig, '-dpng', [figdir, 'marginal_pdf.png']);         
            end
        end
    end
end

function setplotattr(m_plot, chot, clim, eqs_sw, view_angle, shiftorigin)
    colormap(chot);
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
end
% Create colorbar