function params = readParams(filename)

    params = struct();
    fid = fopen(filename,'r');

    if fid == -1
        error('Could not open metadata file: %s', filename);
    end

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Skip empty lines or comments
        if isempty(line) || startsWith(line, {'#','%'})
            continue;
        end

        % Remove inline comments
        cidx = regexp(line, '[#%]', 'once');
        if ~isempty(cidx)
            line = strtrim(line(1:cidx-1));
            if isempty(line)
                continue;
            end
        end

        % Split into key and value
        tokens = split(line, '=');
        key = strtrim(tokens{1});
        val = strtrim(strjoin(tokens(2:end), '='));

        % Convert value
        num = str2double(val);
        if ~isnan(num)
            val = num;
        elseif strcmpi(val,'true')
            val = true;
        elseif strcmpi(val,'false')
            val = false;
        end

        % ---- NEW PART: allow repeated keys ----
        if isfield(params, key)
            if ~iscell(params.(key))
                params.(key) = {params.(key)};
            end
            params.(key){end+1} = val;
        else
            params.(key) = val;
        end
    end

    fclose(fid);

    % Normalize folder to cell array
    if isfield(params,'folder') && ~iscell(params.folder)
        params.folder = {params.folder};
    end
end
