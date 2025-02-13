function [SimData, OrbitData, AttitudeData, estimData] = config_module(name)
    %%                      CONFIGURATION MODULE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This module is in charge of read the configuration
    %   file specified in the input path and load all the parameretes from
    %   orbital, input and attitude.
    %
    %   Input:
    %       name: configuration file name.
    %   
    %   Output: input, orbital and attitude dynamics modules parameters.
    % --------------------------------------------------------------------%

    % add path
    addpath('modules/config_module/functions/');
    
    % path
    path = './config/' + name;

    % Read input variables
    type = "SIMULATION";
    N = 8;
    [SimData] = readConfig(path, type, N);

    % Read orbital variables
    type = "ORBITAL";
    N = 12;
    [OrbitData] = readConfig(path, type, N);

    % Read attitude dynamics varaibles
    type = "ATTITUDE";
    N = 15;
    [AttitudeData] = readConfig(path, type, N);

    % Read estimation parameters
    type = "ESTIMATION";
    N = 8;
    [estimData] = readConfig(path, type, N);
end