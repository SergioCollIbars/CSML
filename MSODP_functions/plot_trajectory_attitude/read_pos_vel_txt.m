function [t, rA, vA, rB, vB]= read_pos_vel_txt(filename)
% read_dual_traj_txt Reads time, position, and velocity for spacecraft A and B.

%

% File format:

%   First line: metadata/header, e.g.

%       2 1201 1202

%

%   Then repeated blocks:

%       time

%       xA yA zA

%       vxA vyA vzA

%       xB yB zB

%       vxB vyB vzB

%

% Outputs:

%   t  : Nx1 time vector

%   rA : Nx3 position matrix for spacecraft A

%   vA : Nx3 velocity matrix for spacecraft A

%   rB : Nx3 position matrix for spacecraft B

%   vB : Nx3 velocity matrix for spacecraft B

    fid = fopen(filename, 'r');

    if fid == -1

        error('Could not open file: %s', filename);

    end

    % Ignore first header line

    fgetl(fid);

    % Read the rest as strings to handle Fortran D exponents

    raw = textscan(fid, '%s');

    fclose(fid);

    raw = raw{1};

    % Convert D exponent to MATLAB-readable E exponent

    raw = strrep(raw, 'D', 'E');

    raw = strrep(raw, 'd', 'E');

    values = str2double(raw);

    % Each block has:

    %   1 time

    %   3 position A

    %   3 velocity A

    %   3 position B

    %   3 velocity B

    % = 13 values

    if mod(numel(values), 13) ~= 0

        error('Unexpected file format: number of values after header is not divisible by 13.');

    end

    n = numel(values) / 13;

    values = reshape(values, 13, n).';

    t  = values(:, 1);

    rA = values(:, 2:4);

    vA = values(:, 5:7);

    rB = values(:, 8:10);

    vB = values(:, 11:13);
end