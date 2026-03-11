function [M] = read_regress(filename)
    fid = fopen(filename,'r');
    assert(fid~=-1, 'Could not open file.');
    
    % 1) Find the "PARTIALS" section
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if strcmp(line,'PARTIALS')
            break
        end
    end
    assert(~feof(fid), 'PARTIALS section not found.');
    
    % 2) Skip description lines until we hit the first numeric row
    %    (and also skip blank lines)
    posDataStart = [];
    while ~feof(fid)
        pos = ftell(fid);
        line = strtrim(fgetl(fid));
        if isempty(line)
            continue
        end
        % If line starts with a number (e.g., "0.0"), it's data
        if ~isempty(regexp(line,'^[\+\-]?\d+(\.\d+)?([eEdD][\+\-]?\d+)?', 'once'))
            posDataStart = pos;
            break
        end
    end
    assert(~isempty(posDataStart), 'Could not find start of numeric PARTIALS data.');
    fseek(fid, posDataStart, 'bof');
    
    % 3) Read numeric rows until "Number of obs" (or non-numeric)
    rows = {};
    while ~feof(fid)
        pos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
    
        if startsWith(s,'Number of obs','IgnoreCase',true)
            break
        end
        if isempty(s)
            continue
        end
    
        % If line is not numeric anymore, stop
        if isempty(regexp(s,'^[\+\-]?\d', 'once'))
            break
        end
    
        % Parse all numbers on the line (handles variable spacing)
        vals = sscanf(s,'%f').';
        rows{end+1,1} = vals; %#ok<SAGROW>
    end
    
    fclose(fid);
    
    % 4) Convert to numeric matrix
    %    (assumes each row has same number of columns)
    nCols = numel(rows{1});
    M = nan(numel(rows), nCols);
    for i = 1:numel(rows)
        M(i,1:numel(rows{i})) = rows{i};
    end
end

