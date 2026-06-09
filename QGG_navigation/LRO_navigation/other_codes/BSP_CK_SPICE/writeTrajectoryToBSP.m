function writeTrajectoryToBSP(t_et, r, v, bspFile, targetID)
% writeTrajectoryToBSP
% Creates a SPICE BSP/SPK file from position and velocity time series.
%
% Inputs:
%   t_et    [1xN] or [Nx1] ET seconds past J2000, TDB
%   r       [3xN] position in km, J2000, relative to Moon
%   v       [3xN] velocity in km/s, J2000, relative to Moon
%   bspFile output BSP filename, e.g. 'my_traj.bsp'

    t_et = t_et(:)';     % make row vector

    if size(r,1) ~= 3
        r = r';
    end

    if size(v,1) ~= 3
        v = v';
    end

    states = [r; v];      % 6xN, km and km/s

    if exist(bspFile, 'file')
        delete(bspFile)
    end

    % SPICE IDs
    centerID = 301;       % Moon

    handle = cspice_spkopn(bspFile, 'MATLAB_TRAJ_SPK', 0);
    step   = t_et(2)-t_et(1);

cspice_spkw08( handle, ...
               targetID, ...
               centerID, ...
               'J2000', ...
               t_et(1), ...
               t_et(end), ...
               'MY_TRAJ', ...
               7, ...
               states, ...
               t_et(1), ...
               step );

    cspice_spkcls(handle);

    fprintf('Created BSP file: %s\n', bspFile);
end