function [outData] = readConfig(path, type, N)
    %%                      readConfig FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    
    %   Date: 26/10/2022
    
    %   Description: Function to read the config files given a path and a
    %   configuration type.
    %
    %   Input:
    %       path: configuration file path.
    %       type: configuration type to read.
    %       N: numbe of variables to read.
    %   
    %   Output: type especified parameters (double) and string messages.
    % --------------------------------------------------------------------%

    % Define output data vector.
    if(type == "ORBITAL")
        outData = zeros(N + 1, 1);
        outData(end) = NaN;
    else
        outData = zeros(N, 1);
    end

    % Point file 
    fileID = fopen(path, 'rt');
    
    % Read module variable
    read = false;
    count = 0;

    % load txt data
    while count < N

        tline = fgetl(fileID);

        if(tline == type)
            read = true;
        elseif(read == true)
             count = count + 1;
             
             a = strsplit(tline, '=');
             b = strsplit(a{2}, ';');
             
             % save data into matrix
             outData(count, 1) = str2double(b{1});
        end
        
    end

    % close file
    fclose(fileID);

end