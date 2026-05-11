function att = read_att_to_rpy(filename)
%READ_ATT_TO_RPY Read MSODP .att file and convert quaternion to roll/pitch/yaw
%
% OUTPUT:
%   att.year   [N x 1]
%   att.doy    [N x 1]
%   att.sod    [N x 1]
%   att.flag1  [N x 1]
%   att.q      [N x 4]   quaternion [q0 q1 q2 q3]
%   att.flag2  [N x 1]
%   att.roll   [N x 1]   rad
%   att.pitch  [N x 1]   rad
%   att.yaw    [N x 1]   rad
%   att.t      [N x 1]   seconds from first record
%
% Assumption:
%   Quaternion columns in file are [q0 q1 q2 q3] (scalar first)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Skip header until end-of-header marker
    line = '';
    while ischar(line)
        line = fgetl(fid);
        if ~ischar(line)
            fclose(fid);
            error('Could not find +eoh______ in file.');
        end
        if contains(line, '+eoh______')
            break;
        end
    end

    % Read numeric records
    % Format from header:
    % (I4,1X,I3,1X,I5,1X,I7,4(1X,E22.15),1X,I8)
    C = textscan(fid, '%f %f %f %f %f %f %f %f %f', ...
        'MultipleDelimsAsOne', true, 'CollectOutput', false);

    fclose(fid);

    year  = C{1};
    doy   = C{2};
    sod   = C{3};
    flag1 = C{4};
% %     q0    = C{5};
% %     q1    = C{6};
% %     q2    = C{7};
% %     q3    = C{8};
    q1    = C{5};
    q2    = C{6};
    q3    = C{7};
    q0    = C{8};

    flag2 = C{9};

    q = [q0, q1, q2, q3];

    % Normalize quaternion just in case
    qnorm = sqrt(sum(q.^2, 2));
    q = q ./ qnorm;

    % Convert to roll-pitch-yaw
    [roll, pitch, yaw] = quat_to_rpy(q);

    % Time from first sample
    t = sod - sod(1) + (doy - doy(1))*86400;

    % Output
    att.year  = year;
    att.doy   = doy;
    att.sod   = sod;
    att.flag1 = flag1;
    att.q     = q;
    att.flag2 = flag2;
    att.roll  = roll;
    att.pitch = pitch;
    att.yaw   = yaw;
    att.t     = t;
end