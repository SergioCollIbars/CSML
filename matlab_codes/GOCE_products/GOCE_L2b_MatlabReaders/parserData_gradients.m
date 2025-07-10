clear;
clc;
close all;

folderPath = 'input_data';  % ← Replace with your folder

% Get list of all .EDF files in the folder
edfFiles = dir(fullfile(folderPath, '*.EDF'));

% output data
dataBlock = ones(86400*40, 7) * NaN;
posVec = 1;
for i = 1:length(edfFiles)
    fileName = edfFiles(i).name;

    if ~contains(fileName, 'EGG_NOM')
        continue;
    end

    filePath = fullfile(folderPath, fileName);
    fid = fopen(filePath, 'r');
    disp("Reading EGG_NOM: " + string(filePath));

    count = 1;
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line)
            continue;
        end

        if count == 1
            count = count + 1;
            continue;
        end

        % Read full line of floats
        dataRow = sscanf(line, '%f')';
        if length(dataRow) < 5
            continue; % skip malformed rows
        end

        gpsTime = dataRow(1);
        gradients = dataRow(2:7); % xx, yy, zz, xy, xz, yz
        dataBlock(posVec, :) = [gpsTime, gradients];
        
        posVec = posVec + 1;
    end

    fclose(fid);
end

% plot gradients
idx = find(~isnan(dataBlock(:, 1)));
data = dataBlock(idx, :);

t = data(:, 1);
gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
t_UTC = gps_epoch + seconds(t);        % date time 


figure()
tt  = ["\Gamma_{xx}", "\Gamma_{yy}", "\Gamma_{zz}", ...
       "\Gamma_{xy}", "\Gamma_{xz}", "\Gamma_{yz}"];
for j = 1:6
    subplot(2, 3, j);
    plot(t_UTC, data(:, j+1)./1E-9, 'LineWidth', 2);
    ylabel('Eotvos');
    title(tt(j));
    grid on;
end