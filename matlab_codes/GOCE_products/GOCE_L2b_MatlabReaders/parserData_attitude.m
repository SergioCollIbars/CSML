clear;
clc;
close all;

folderPath = 'input_data';
edfFiles = dir(fullfile(folderPath, '*.EDF'));

% GPS Epoch (1980-01-06 00:00:00 UTC)
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);

% output data
dataBlock = [];

for i = 1:length(edfFiles)
    filePath = fullfile(folderPath, edfFiles(i).name);
    fid = fopen(filePath, 'r');
    disp('Reading ...' + string(filePath));

    firstEpochStr = '';
    readingHeader = true;
    
    while ~feof(fid)
        line = fgetl(fid);
        
        % Skip empty lines
        if isempty(line)
            continue;
        end
        
        % Extract the first epoch if found
        if startsWith(line, '# First epoch:')
            firstEpochStr = strtrim(erase(line, '# First epoch:'));
            dt = datetime(firstEpochStr);
            gpsSeconds = seconds(dt - gpsEpoch);
        end
        
        % Check for end of header
        if startsWith(line, '# End of header')
            readingHeader = false;
            continue;
        end
        
        % Read data rows after header ends
        if ~readingHeader && ~startsWith(line, '#')
            dataRow = sscanf(line, '%f')';  % Transpose to row vector
            dataRow(1) = dataRow(1) + gpsSeconds;
            dataBlock = [dataBlock; dataRow]; %#ok<AGROW>
        end
    end
    
    fclose(fid);
end

matFilePath = fullfile("data/", 'Nov_L2position.mat');

if exist(matFilePath, 'file') == 2
    % File exists — perform your action
    disp('MAT file exists. Loading data...');
    load(matFilePath);

    % get common GPS times
    [commonTimes, idx1, idx2] = intersect(positions(:, 1), dataBlock(:, 1));

    % output file
    outPut = dataBlock(idx2, :);

else
    % File does not exist — perform a different action
    disp('MAT file not found. Processing raw EDF files...');
    
    
    % output file
    outPut = dataBlock;
end

% save data into arrays 
name ="Nov_L2ECEF2ITRF.mat";
save("data/"+name, "outPut");
