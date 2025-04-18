clear;
clc;
close all;

filename = 'GO_CONS_SST_PKI_2__20121115T235944_20121116T235943_0201.EDF';

fid = fopen(filename);

% Initialize output arrays
positions = [];

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
    end
end

fclose(fid);

% save data into arrays 
dateStr = datestr(dt, 'dd_mmm_yyyy');
name = dateStr + "_L2position.mat";
save("data/"+name, "positions");