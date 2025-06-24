clear;
clc;
close all;

disp('Read files ...')
folderPath = 'input_data';
edfFiles = dir(fullfile(folderPath, '*.EDF'));

% GPS Epoch (1980-01-06 00:00:00 UTC)
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);

% output data
dataBlock = ones(86400*40, 5) * NaN;

count = 1;
for i = 1:length(edfFiles)
    fileName = edfFiles(i).name;
    
    % Only process files that contain 'SST_PRM'
    if ~contains(fileName, 'SST_PRM')
        continue;
    end

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
            % % dataBlock = [dataBlock; dataRow]; %#ok<AGROW>
            dataBlock(count, :) = dataRow;
            count = count + 1;
        end
    end
    
    fclose(fid);
end

% remove extra vector rows
idx = ~isnan(dataBlock(:, 1));
dataBlock = dataBlock(idx, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2Block = ones(86400*40, 5) * NaN;
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
        lastFour = dataRow(end-3:end);
        data2Block(posVec, :) = [gpsTime, lastFour];
        
        posVec = posVec + 1;
    end

    fclose(fid);
end

% remove extra vector rows
idx = ~isnan(data2Block(:, 1));
data2Block = data2Block(idx, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Interpolate quaternions ....')

if(~isempty(data2Block))
    % match common GPS times (integreer)
    time_block1 = dataBlock(:, 1); time_block2 = floor(data2Block(:, 1));
    [~, idx1, idx2] = intersect(time_block1, time_block2);
    
    dataBlock = dataBlock(idx1, :); data2Block = data2Block(idx2, :);
    
    % interpolate quaternions to the data block 1 times.
    for j = 2:length(dataBlock(:, 1))
         qa = data2Block(j-1, 2:end); ta = data2Block(j-1, 1);
         qb = data2Block(j, 2:end);   tb = data2Block(j, 1);
         t = dataBlock(j, 1);
         qt = interpolateQuaternionGOCE(qa, qb, ta, tb, t);
         data2Block(j, 2:end) = qt;
         data2Block(j, 1)     = t;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save Data ...')

matFilePath = fullfile("data/", 'Nov_L2position.mat');

if exist(matFilePath, 'file') == 2
    % File exists — perform your action
    disp('MAT file exists. Loading data...');
    load(matFilePath);
    
    t1 = positions(:, 1); t2  = dataBlock(:, 1);
    t3 = data2Block(:, 1);
    commonTimes = intersect(intersect(t1, t2), t3);

    [~, idx2] = ismember(commonTimes, t2);
    [~, idx3] = ismember(commonTimes, t3);

    % get common GPS times
    outPut  = dataBlock(idx2, :);
    outPut2 = data2Block(idx3, :);

else
    % File does not exist — perform a different action
    disp('MAT file not found. Processing raw EDF files...');
    
    
    % output file
    outPut  = dataBlock;
    outPut2 = data2Block;
 end

% save data into arrays 
name ="Nov_L2ECEF2ITRF.mat";
save("data/"+name, "outPut");

name ="Nov_L2ITRF2GRF.mat";
save("data/"+name, "outPut2");


%% FUNCTIONS
function qt = interpolateQuaternionGOCE(qa, qb, ta, tb, t)
%INTERPOLATEQUATERNIONGOCE Interpolates between two unit quaternions qa and qb at time t
% using the method from the GOCE Level 2 Product Data Handbook (Equations 4.10–4.17).
%
% Inputs:
%   qa : 1x4 unit quaternion at time ta (format: [q1 q2 q3 q4])
%   qb : 1x4 unit quaternion at time tb (format: [q1 q2 q3 q4])
%   ta : scalar, time of qa
%   tb : scalar, time of qb
%   t  : scalar, desired interpolation time (ta <= t <= tb)
%
% Output:
%   qt : 1x4 interpolated quaternion at time t (unit quaternion)

% Ensure row vectors
qa = qa(:)';
qb = qb(:)';

% Step 1: Sign correction (Eq. 4.10)
dotVec = dot(qa(1:3), qb(1:3));
if dotVec < 0
    qb = -qb;
end

% Step 2: Compute differential quaternion q_ab = qa^* * qb (Eq. 4.11)
% qa* = [ -qa(1:3), qa(4) ]
qa_conj = [-qa(1:3), qa(4)];
q_ab = quatMultiply(qa_conj, qb);

% Step 3: Rotation angle Φ_ab = 2 * acos(q_ab4) (Eq. 4.13)
phi_ab = 2 * acos(q_ab(4));

% Interpolation not needed if angle is zero (qa == qb)
if phi_ab < 1e-10
    qt = qa;
    return;
end

% Step 4: Linearly interpolate angle Φ_at (Eq. 4.14)
phi_at = phi_ab * (t - ta) / (tb - ta);

% Step 5: Compute interpolated rotation quaternion q_at (Eq. 4.15)
q_at = zeros(1,4);
q_at(4) = cos(phi_at / 2);
scale = sin(phi_at / 2) / sin(phi_ab / 2);
q_at(1:3) = q_ab(1:3) * scale;

% Step 6: Final interpolated quaternion qt = qa * q_at (Eq. 4.16–4.17)
qt = quatMultiply(qa, q_at);

% Normalize to avoid numerical drift
qt = qt / norm(qt);
end

% Helper function: Quaternion multiplication
function q = quatMultiply(q1, q2)
% Multiply two quaternions: q = q1 * q2
% Each q is a 1x4 vector: [q1 q2 q3 q4] (q4 is scalar part)

x1 = q1(1); y1 = q1(2); z1 = q1(3); w1 = q1(4);
x2 = q2(1); y2 = q2(2); z2 = q2(3); w2 = q2(4);

x =  w1*x2 + x1*w2 + y1*z2 - z1*y2;
y =  w1*y2 - x1*z2 + y1*w2 + z1*x2;
z =  w1*z2 + x1*y2 - y1*x2 + z1*w2;
w =  w1*w2 - x1*x2 - y1*y2 - z1*z2;

q = [x, y, z, w];
end



