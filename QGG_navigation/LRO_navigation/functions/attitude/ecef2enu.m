function ENU_ECEF = ecef2enu(r0)
% r0: 3x1 ECEF position of ENU origin [m]
x = r0(1); y = r0(2); z = r0(3);

lon = atan2(y, x);
lat = atan2(z, hypot(x,y));   % geocentric latitude

sl = sin(lon);  cl = cos(lon);
sp = sin(lat);  cp = cos(lat);

ENU_ECEF = [ -sl,        cl,      0;
      -sp*cl, -sp*sl,   cp;
       cp*cl,  cp*sl,   sp ];
end