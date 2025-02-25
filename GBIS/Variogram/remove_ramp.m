function [] = remove_ramp(inputFile, wavelength, pos_path)
    close all
    clc
    if nargin < 3
        pos_path = [];
    end

    % Load dataset
    disp('Ingesting data to estimate (semi-)variogram ...')
    insarData = load(inputFile);

    % Create colormap
    cmapSeismo = colormap_cpt('GMT_seis.cpt', 256);    % GMT 'Seismo' colormap for wrapped data

    % Find a local reference point
    refPoint = [min(insarData.Lon), min(insarData.Lat)]; % Determine local reference point

    % Convert phase to LOS displacement
    convertedPhase = (insarData.Phase / (4 * pi)) * wavelength;  % Convert phase from radians to m
    los = single(-convertedPhase);                              % Convert to Line-of-sight displacement in m

    % Determine subsampling factor for faster plotting
    if length(los) > 400000 && length(los) < 1000000
        sampling = 2;
    elseif length(los) > 1000000
        sampling = 5;
    else
        sampling = 1;
    end

    % Extract subset from rectangular area or after masking
    subset = los;
    llon = insarData.Lon;
    llat = insarData.Lat;
    inc = insarData.Inc;
    heading = insarData.Heading;

    if ~isempty(pos_path)
        % Use the provided polygon coordinates
        load(pos_path);
        if exist('pos', 'var') && ~isempty(pos)
            disp('Using provided polygon coordinates.');
            in = inpolygon(insarData.Lon, insarData.Lat, pos(:,1), pos(:,2));
            ixSubset = find(in == 0); % Points outside the polygon
        elseif exist('pos_rectangular', 'var') && ~isempty(pos_rectangular)
            disp('Using provided rectangular coordinates.');
            % Extract the min and max values from the saved rectangular coordinates
            maxLon = pos_rectangular(1);
            maxLat = pos_rectangular(2);
            minLon = pos_rectangular(3);
            minLat = pos_rectangular(4);
            ixSubset = find(insarData.Lat > minLat & insarData.Lat < maxLat & insarData.Lon < minLon & insarData.Lon > maxLon);
        else
            error('The provided file does not contain valid polygon or rectangular coordinates.');
        end
    end

    % Subset points outside the polygon/rectangle
    subset_outside = los(ixSubset);
    llon_outside = insarData.Lon(ixSubset);
    llat_outside = insarData.Lat(ixSubset);
    inc_outside = insarData.Inc(ixSubset);
    heading_outside = insarData.Heading(ixSubset);

    % Subset points inside the polygon/rectangle
    inside_indices = ~ismember(1:length(los), ixSubset);
    subset_inside = los(inside_indices);  % Complement of ixSubset
    llon_inside = insarData.Lon(inside_indices);
    llat_inside = insarData.Lat(inside_indices);
    inc_inside = insarData.Inc(inside_indices);
    heading_inside = insarData.Heading(inside_indices);

    % Apply the 3-sigma rule to subset_outside
    mu = nanmean(subset_outside);
    sigma = nanstd(subset_outside);
    outliers = subset_outside > (mu + 3 * sigma) | subset_outside < (mu - 3 * sigma);
    % Remove outliers from subset_outside
    subset_outside = subset_outside(~outliers);
    llon_outside = llon_outside(~outliers);
    llat_outside = llat_outside(~outliers);
    inc_outside = inc_outside(~outliers);
    heading_outside = heading_outside(~outliers);

    % Merge inside and cleaned outside data for all fields
    merged_subset = [subset_inside; subset_outside];
    merged_llon = [llon_inside; llon_outside];
    merged_llat = [llat_inside; llat_outside];
    merged_inc = [inc_inside; inc_outside];
    merged_heading = [heading_inside; heading_outside];

    % Display subregion (non-detrended)
    figure('Position', [1, 1, 1200, 1000]);
    subplot(1, 3, 1)
    scatter(merged_llon(:), merged_llat(:), [], merged_subset(:), '.');
    colormap(cmapSeismo);
    caxis([-1 1]);
    axis xy
    axis equal
    axis tight
    title('NON-DETRENDED')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    colorbar

    %% Remove linear trend from the merged subregion
    sll = [merged_llon'; merged_llat']';
    xy = llh2local(sll', refPoint);
    xy = xy * 1000;

    A = [xy' ones([length(xy) 1])];   % Design matrix for linear ramp removal

    % Generate the corresponding .cov filename
    [path, name, ~] = fileparts(inputFile);   % Extract base name of input file (without extension)
    covFileName = fullfile(path, [name, '.cov']);  % Append .cov to the base name

    % Load the coefficients from the .cov file
    fileID = fopen(covFileName, 'r');
    if fileID == -1
        error('Cannot open the .cov file: %s', covFileName);
    end

    % Read through the .cov file to extract the last three coefficients
    coeffs = [];
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, 'Coeff')
            coeffValue = sscanf(line, 'Coeff: %f');
            coeffs = [coeffs; coeffValue];
        end
    end
    fclose(fileID);

    % Ensure coeffs contains the expected three values
    if length(coeffs) ~= 3
        error('Loaded coefficient data from .cov does not contain 3 values.');
    end

    coeff = coeffs;
    deramped = merged_subset - A * coeff;  % Remove the linear trend

    % Convert the deramped LOS displacement back to Phase (in radians)
    derampedPhase = (-deramped / wavelength) * (4 * pi);

    % Replace the updated data in insarData with the new merged values
    insarData.Lon = merged_llon;
    insarData.Lat = merged_llat;
    insarData.Phase = derampedPhase;  % Update Phase with detrended values
    insarData.Inc = merged_inc;  % Update Incidence
    insarData.Heading = merged_heading;  % Update Heading

    % Save the new .mat file with the updated data
    saveFileName = fullfile(path, [name, '_rr.mat']);
    save(saveFileName, '-struct', 'insarData');
    disp(['Removed ramp values saved to: ', saveFileName]);

    %% Display the trend and subregion after trend removal
    subplot(1, 3, 2)
    scatter(merged_llon(:), merged_llat(:), [], A(:, :) * coeff, '.');
    colormap(cmapSeismo);
    caxis([-0.5 0.5]);
    axis xy
    axis equal
    axis tight
    title('ESTIMATED TREND')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    colorbar

    subplot(1, 3, 3)
    scatter(merged_llon(:), merged_llat(:), [], deramped(:), '.');
    colormap(cmapSeismo);
    caxis([-1 1]);
    axis xy
    axis equal
    axis tight
    title('DETRENDED')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    colorbar

    % Save the figure as a PNG file
    pngFileName = fullfile(path, [name, '.removed.ramp.png']);
    saveas(gcf, pngFileName);
    disp(['Figure saved as: ', pngFileName]);

end
