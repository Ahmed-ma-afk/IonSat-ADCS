function [latitude, longitude, altitude] = ecef2geod(x, y, z)

% ECEF2GEOD Convert ECEF coordinates to geodetic coordinates.
% 
% Usage: [LATITUDE, LONGITUDE, ALTITUDE] = ECEF2GEOD(X, Y, Z, TOL)
% 
% Converts Earth-centered, Earth fixed (ECEF) coordinates X, Y, and Z to
% geodetic coordinates LATITUDE, LONGITUDE, and ALTITUDE. For a matrix
% input, the first dimension with length 3 is assumed to have the three
% separate X, Y, and Z inputs across it. The World Geodetic System 1984
% (WGS84) ellipsoid model of the Earth is assumed.
% 
% Inputs:
%   -X: x coordinates of the point in meters.
%   -Y: y coordinates of the point in meters.
%   -Z: z coordinates of the point in meters.
% 
% Ouputs:
%   -LATITUDE: Geodetic latitude in degrees.
%   -LONGITUDE: Geodetic longitude in degrees.
%   -ALTITUDE: Height above the Earth in meters.

tol=1e-12;
% WGS84 parameters.
a = 6378137; f = 1/298.257223563; b = a*(1 - f); e2 = 1 - (b/a)^2;

% Longitude is easy:
longitude = atan2(y, x)*180/pi;

% Compute latitude recursively.
rd = hypot(x, y);
[latitude, Nphi] = recur(asin(z ./ hypot(x, hypot(y, z))), z, a, e2, ...
    rd, tol, 1);
sinlat = sin(latitude); coslat = cos(latitude); 
latitude = latitude*180/pi;

% Get altitude from latitude.
altitude = rd.*coslat + (z + e2*Nphi.*sinlat).*sinlat - Nphi;


function [latitude, Nphi] = recur(lat_in, z, a, e2, rd, tol, iter)
thisNphi = a ./ sqrt(1 - e2*sin(lat_in).^2);
nextlat = atan((z + thisNphi*e2.*sin(lat_in))./rd);
if all(abs(lat_in - nextlat) < tol) || iter > 100
    latitude = nextlat; Nphi = thisNphi;
else
    [latitude, Nphi] = recur(nextlat, z, a, e2, rd, tol, iter + 1);
end