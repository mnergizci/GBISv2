% Define the directory containing the .mat files
inputDir = pwd;  % Replace with the actual path if needed

% Get a list of all .mat files in the directory
matFiles = dir(fullfile(inputDir, '*.mat'));

% Loop through each file, filter based on 'azi' or 'boi', and process them
for i = 1:length(matFiles)
    % Get the file name
    fileName = matFiles(i).name;
    
    % Check if the file name contains 'azi' or 'boi'
    if contains(fileName, 'azi', 'IgnoreCase', true) || contains(fileName, 'boi', 'IgnoreCase', true)
        % Load the .mat file
        filePath = fullfile(inputDir, fileName);
        data = load(filePath);
        
        % Check if the file contains 'Heading' variable
        if isfield(data, 'Heading')
            % Check the 4th character of the file name to determine adjustment
            if length(fileName) >= 4
                fourthChar = fileName(4);

                % Adjust the Heading based on the 4th character
                if fourthChar == 'A'
                    % Subtract 180 from Heading if the 4th character is 'A'
                    data.Heading = data.Heading - 180;
                    disp(['Subtracted 180 from Heading for file: ' fileName]);
                elseif fourthChar == 'D'
                    % Add 180 to Heading if the 4th character is 'D'
                    data.Heading = data.Heading + 180;
                    disp(['Added 180 to Heading for file: ' fileName]);
                else
                    % Display a message if the 4th character is not 'A' or 'D'
                    disp(['No adjustment needed for file: ' fileName ' (4th character is not A or D)']);
                end
                
                % Save the updated variable back to the .mat file
                save(filePath, '-struct', 'data');
            else
                % Display a message if the file name is too short
                disp(['File name is too short to check 4th character: ' fileName]);
            end
        else
            % Display a message if the 'Heading' variable is not found
            disp(['Heading variable not found in file: ' fileName]);
        end
    end
end
