function [ECI_ECEF] = convert_ECEF2ECI(gps_seconds)
    % convert to Julian Date (JD)
    gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
    utc_time = gps_epoch + seconds(gps_seconds);
    JD = juliandate(utc_time);

    % Time in Julian centuries from J2000
    T = (JD - 2451545.0)/36525;
    
    % GMST in seconds
    GMST_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T + ...
               0.093104*T^2 - 6.2e-6*T^3;
    
    % Convert to radians and wrap to [0, 2pi]
    GMST = mod(deg2rad(mod(GMST_sec/240, 360)), 2*pi); % divide by 240 to convert seconds to degrees

    % rotation matrix
    ECI_ECEF = [ cos(GMST), sin(GMST), 0;
                -sin(GMST), cos(GMST), 0;
                  0,         0, 1];

end

