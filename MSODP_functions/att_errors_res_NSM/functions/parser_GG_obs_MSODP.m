function [out] = parser_GG_obs_MSODP(folderPath)
    %PARSER_GG_OBS_MSODP Summary of this function goes here
    %   Detailed explanation goes here
    
    % Check folder exists
    if ~isfolder(folderPath)
        error('Folder does not exist.');
    end

    % Get list of all .txt files
    files = dir(fullfile(folderPath, '*.ggr'));

    if isempty(files)
        warning('No .txt files found in folder.');
        return;
    end

    pat = '(?<date>\d{4}-\d{2}-\d{2})';
    t0 = datetime(2000,1,1,12,0,0,'TimeZone','UTC'); % "J2000" in UTC approx
    
    dataOutput = nan(60*86400, 7);
    firstIndx  = 1;
    for k = 1:length(files)

        fileName = files(k).name;
        fullPath = fullfile(folderPath, fileName);

        tok = regexp(fileName, pat, 'names', 'once');
        dstr= string(tok.date);
        t   = datetime(char(dstr),'InputFormat','yyyy-MM-dd','TimeZone','UTC');
        secJ2K = seconds(t - t0);

        fprintf('Loading: %s\n', fileName);
         fid = fopen(fullPath, 'r');
        
         if fid == -1
            error('Cannot open file %s', fileName);
        end

        C = textscan(fid, '%f', ...
                     'Delimiter', {' ', '\t'}, ...
                     'MultipleDelimsAsOne', true);

        fclose(fid);

        % Convert long column vector into matrix
        raw = C{1};

        % Determine number of columns from first line
        firstLine = fileread(fullPath);
        firstLine = extractBefore(firstLine, newline);
        nCols = numel(strsplit(strtrim(firstLine)));

        data          = reshape(raw, nCols, []).';
        time_day_sec  = data(:, 3) + secJ2K;
        lastIndx      =  length(data) + firstIndx -1;
        dataOutput(firstIndx:lastIndx, :) = [time_day_sec, data(:, 21:26)];

        firstIndx = lastIndx + 1;
    end
    
    % remove nan index
    idx = find(~isnan(dataOutput(:, 1)));
    out = dataOutput(idx, :); 
end

