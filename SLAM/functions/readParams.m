function params = readParams(filename)
    params = struct();
    fid = fopen(filename,'r');

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Skip empty lines or comments
        if isempty(line) || startsWith(line, {'#','%'})
            continue;
        end

        % Split into key and value
        tokens = split(line, '=');
        key = strtrim(tokens{1});
        val = strtrim(strjoin(tokens(2:end), '='));

        % Convert numeric values if possible
        num = str2double(val);
        if ~isnan(num)
            params.(key) = num;
        elseif strcmpi(val,'true')
            params.(key) = true;
        elseif strcmpi(val,'false')
            params.(key) = false;
        else
            params.(key) = val;  % String
        end
    end

    fclose(fid);
end
