function fault_patches_generator(filename, orig, mode)
    % This code is improved to divide the fault segments into variable
    % patches:
    %   if the mode equals '1' the patch sizes = 1-3-6-10 km along the downdip.
    % MNergizci, 2024

    % Open the file
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('File cannot be opened: %s', filename);
    end

    % Read the first line (header)
    header = fgetl(fileID);

    % Initialize segment count
    fault_segment = 0;

    % Create an empty cell array to store each line of data
    data = {};

    % Read the file line by line
    while ~feof(fileID)
        line = fgetl(fileID);
        % Check if the line is not empty
        if ischar(line) && ~isempty(line)
            fault_segment = fault_segment + 1;
            data{fault_segment} = line; % Store the line of data
        end
    end

    % Close the file
    fclose(fileID);

    % Display the number of segments
    disp(' ');
    disp(['Number of segments: ', num2str(fault_segment)]);
    disp(' ');

    % Initialize the structure array
    patch(1:fault_segment) = struct('lat', 0, 'lon', 0, 'dip', 0, 'strike', 0, 'length', 0, 'depth', 0, 'size_patches', 0, 'inc', 0, 'X', 0, 'Y', 0, 'Z', 0);

    % Parse the data and store in the structure array
    for i = 1:fault_segment
        % Split the line by commas
        values = str2num(data{i});

        % Assign the values to the corresponding fields in the structure
        patch(i).lat = values(1);
        patch(i).lon = values(2);
        patch(i).dip = values(3);
        patch(i).strike = values(4);
        patch(i).length = values(5) * 1000; % Convert km to meters
        patch(i).depth = values(6) * 1000;  % Convert km to meters
        patch(i).size_patches = values(7) * 1000; % Convert km to meters
        patch(i).inc = values(8);
    end

    % Define the origin coordinates
    patch_ll = [patch.lon; patch.lat];
    patch_XY = llh2local([patch_ll; zeros(1, fault_segment)], orig) * 1000;

    % Assign the X and Y values from patch_XY to the corresponding fields in the structure
    for i = 1:fault_segment
        patch(i).X = patch_XY(1, i);
        patch(i).Y = patch_XY(2, i);
        patch(i).Z = zeros(1); % Placeholder for Z field
    end
    
    % Display dip, strike, total length, and total width, and patch sizes if mode is 1
    if mode == 1
        m= [];
    
        for i = 1:fault_segment
            dip = patch(i).dip * pi / 180; % Convert to radians
            strike = patch(i).strike * pi / 180; % Convert to radians
            total_length = patch(i).length;
            total_width = patch(i).depth/sin(-1*dip);
            X_cent = patch(i).X;
            Y_cent = patch(i).Y;


            init_patch_sizes = [1, 3, 6, 10]; % Initial patch sizes in km
            init_patch_with_dip=init_patch_sizes./sin(-1*dip)
            cumulative_sum = cumsum(init_patch_with_dip); % [1, 4, 10, 20]
            cum_top_patch = [0, cumulative_sum(1:end-1)]; % [0, 1, 4, 10]
    
            disp(' ');
            disp(['Fault Segment-', num2str(i)]);
            disp(['Dip: ', num2str(dip)]);
            disp(['Strike: ', num2str(strike)]);
            disp(['Total Length: ', num2str(total_length)]);
            disp(['Total Width: ', num2str(total_width)]);
                
            % Compute patch sizes
            for j = 1:length(init_patch_sizes)
                ps = init_patch_sizes(j);
                ps_m = ps * 1000; % Convert to meters
                ps_m_dip=ps_m/sin(-1*dip);
                patch_num_strike = round(total_length / ps_m);
                patch_size_strike = total_length / patch_num_strike;

    
                disp(['Patch size along strike: ', num2str(patch_size_strike), ';   Patch size along dip: ', num2str(ps_m_dip)]);
                    
                % Create the meshgrid
                [AS, ~] = meshgrid([0:patch_size_strike:(patch_num_strike-1)*patch_size_strike] - (patch_num_strike-1)*patch_size_strike/2, [0]);
                UD = ones(size(AS)) * cum_top_patch(j) * 1000; % Fill UD with the constant value
                
                m_out = zeros(10, patch_num_strike);
                m_out(1,:) = repmat(patch_size_strike,1,patch_num_strike);
                m_out(2,:) = repmat(ps_m_dip,1,patch_num_strike);
                m_out(3,:) = repmat(0,1,patch_num_strike) - UD*sin(dip);
                m_out(4,:) = repmat(dip*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(5,:) = repmat(strike*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(6,:) = AS*sin(strike) - UD*cos(strike)*cos(dip) + X_cent;
                m_out(7,:) = AS*cos(strike) + UD*sin(strike)*cos(dip) + Y_cent;
                m = [m, m_out];
            end
        end

        % Save the m_total variable as a .mat file
        save('m_total_mode1.mat', 'm');
        
        % Plotting
        figure
        drawmodel(m)
        axis;
    end
    
    
        % Display dip, strike, total length, and total width, and patch sizes if mode is 1
    if mode == 2
        m= [];
    
        for i = 1:fault_segment
            dip = patch(i).dip * pi / 180; % Convert to radians
            strike = patch(i).strike * pi / 180; % Convert to radians
            total_length = patch(i).length;
            total_width = patch(i).depth/sin(-1*dip);
            X_cent = patch(i).X;
            Y_cent = patch(i).Y;


            init_patch_sizes = [2, 4, 6, 8]; % Initial patch sizes in km
            init_patch_with_dip=init_patch_sizes./sin(-1*dip)
            cumulative_sum = cumsum(init_patch_with_dip); % [1, 4, 10, 20]
            cum_top_patch = [0, cumulative_sum(1:end-1)]; % [0, 1, 4, 10]
    
            disp(' ');
            disp(['Fault Segment-', num2str(i)]);
            disp(['Dip: ', num2str(dip)]);
            disp(['Strike: ', num2str(strike)]);
            disp(['Total Length: ', num2str(total_length)]);
            disp(['Total Width: ', num2str(total_width)]);
                
            % Compute patch sizes
            for j = 1:length(init_patch_sizes)
                ps = init_patch_sizes(j);
                ps_m = ps * 1000; % Convert to meters
                ps_m_dip=ps_m/sin(-1*dip);
                patch_num_strike = round(total_length / ps_m);
                patch_size_strike = total_length / patch_num_strike;

    
                disp(['Patch size along strike: ', num2str(patch_size_strike), ';   Patch size along dip: ', num2str(ps_m_dip)]);
                    
                % Create the meshgrid
                [AS, ~] = meshgrid([0:patch_size_strike:(patch_num_strike-1)*patch_size_strike] - (patch_num_strike-1)*patch_size_strike/2, [0]);
                UD = ones(size(AS)) * cum_top_patch(j) * 1000; % Fill UD with the constant value
                
                m_out = zeros(10, patch_num_strike);
                m_out(1,:) = repmat(patch_size_strike,1,patch_num_strike);
                m_out(2,:) = repmat(ps_m_dip,1,patch_num_strike);
                m_out(3,:) = repmat(0,1,patch_num_strike) - UD*sin(dip);
                m_out(4,:) = repmat(dip*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(5,:) = repmat(strike*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(6,:) = AS*sin(strike) - UD*cos(strike)*cos(dip) + X_cent;
                m_out(7,:) = AS*cos(strike) + UD*sin(strike)*cos(dip) + Y_cent;
                m = [m, m_out];
            end
        end

        % Save the m_total variable as a .mat file
        save('m_total_mode2.mat', 'm');
        
        % Plotting
        figure
        drawmodel(m)
        axis;
    end
    if mode == 3
        m= [];
    
        for i = 1:fault_segment
            dip = patch(i).dip * pi / 180; % Convert to radians
            strike = patch(i).strike * pi / 180; % Convert to radians
            total_length = patch(i).length;
            total_width = patch(i).depth/sin(-1*dip);
            X_cent = patch(i).X;
            Y_cent = patch(i).Y;


            init_patch_sizes = [2, 4, 6, 8, 10]; % Initial patch sizes in km
            init_patch_with_dip=init_patch_sizes./sin(-1*dip)
            cumulative_sum = cumsum(init_patch_with_dip); % [1, 4, 10, 20]
            cum_top_patch = [0, cumulative_sum(1:end-1)]; % [0, 1, 4, 10]
    
            disp(' ');
            disp(['Fault Segment-', num2str(i)]);
            disp(['Dip: ', num2str(dip)]);
            disp(['Strike: ', num2str(strike)]);
            disp(['Total Length: ', num2str(total_length)]);
            disp(['Total Width: ', num2str(total_width)]);
                
            % Compute patch sizes
            for j = 1:length(init_patch_sizes)
                ps = init_patch_sizes(j);
                ps_m = ps * 1000; % Convert to meters
                ps_m_dip=ps_m/sin(-1*dip);
                patch_num_strike = round(total_length / ps_m);
                patch_size_strike = total_length / patch_num_strike;

    
                disp(['Patch size along strike: ', num2str(patch_size_strike), ';   Patch size along dip: ', num2str(ps_m_dip)]);
                    
                % Create the meshgrid
                [AS, ~] = meshgrid([0:patch_size_strike:(patch_num_strike-1)*patch_size_strike] - (patch_num_strike-1)*patch_size_strike/2, [0]);
                UD = ones(size(AS)) * cum_top_patch(j) * 1000; % Fill UD with the constant value
                
                m_out = zeros(10, patch_num_strike);
                m_out(1,:) = repmat(patch_size_strike,1,patch_num_strike);
                m_out(2,:) = repmat(ps_m_dip,1,patch_num_strike);
                m_out(3,:) = repmat(0,1,patch_num_strike) - UD*sin(dip);
                m_out(4,:) = repmat(dip*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(5,:) = repmat(strike*180/pi,1,patch_num_strike); % Convert back to degrees
                m_out(6,:) = AS*sin(strike) - UD*cos(strike)*cos(dip) + X_cent;
                m_out(7,:) = AS*cos(strike) + UD*sin(strike)*cos(dip) + Y_cent;
                m = [m, m_out];
            end
        end

        % Save the m_total variable as a .mat file
        save('m_total_mode3.mat', 'm');
        
        % Plotting
        figure
        drawmodel(m)
        axis;
    end
end
