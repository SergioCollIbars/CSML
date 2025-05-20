clear;
clc;
close all;

folderPath = 'input_data';  % ← Replace with your folder

% Get list of all .EDF files in the folder
edfFiles = dir(fullfile(folderPath, '*.EDF'));

% Initialize output arrays
positions = []; velocity = [];

for k = 1:length(edfFiles)
    filename = fullfile(folderPath, edfFiles(k).name);
    fid = fopen(filename);
    disp('Processing ... ' + string(filename));

    if fid == -1
        warning('Failed to open file: %s', filename);
        continue;
    end

    while ~feof(fid)
        line = fgetl(fid);
        
        if startsWith(line, '*')  % Epoch line
            data = sscanf(line, '* %d %d %d %d %d %f');
            year    = data(1);
            month   = data(2);
            day     = data(3);
            hour    = data(4);
            minute  = data(5);
            second  = data(6);
    
            dt = datetime(year, month, day, hour, minute, second);
            
            % GPS Epoch (1980-01-06 00:00:00 UTC)
            gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
    
            % Calculate the difference in seconds from GPS epoch
            gpsSeconds = seconds(dt - gpsEpoch);
            
        elseif startsWith(line, 'P')  % Position line
            % Satellite PRN = line(2:4); % optional
            tokens = regexp(line, '[\+\-]\d+\.\d+', 'match');
            numbers = str2double(tokens);
            X = numbers(1) * 1000; % km to meters
            Y = numbers(2) * 1000;
            Z = numbers(3) * 1000;
            
            % Save the data
            positions(end+1,:) = [gpsSeconds X Y Z];

        elseif startsWith(line, 'V')  % Velocity line
             % Satellite PRN = line(2:4); % optional
            tokens = regexp(line, '[\+\-]\d+\.\d+', 'match');
            numbers = str2double(tokens);
            X = numbers(1) * 0.1; % decimeter /sec to meters /sec
            Y = numbers(2) * 0.1;
            Z = numbers(3) * 0.1;
            
            % Save the data
            velocity(end+1,:) = [gpsSeconds X Y Z];
        end
    end

    fclose(fid);
end


% save data into arrays 
name ="Nov_L2position.mat";
save("data/"+name, "positions");

name ="Nov_L2velocity.mat";
save("data/"+name, "velocity");