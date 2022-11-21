function [B_ECEF] = NED2ecef(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Convert from spherical to ECEF


% % Convert from spherical to (x,y,z) = (North,East,Down).
% Bx = -Bt;
% By = Bp;
% Bz = -Br;
Bt=-u(1);
Bp=u(2);
Br=-u(3);
lat=u(4);
lon=u(5);

phi = deg2rad(lon); %phi is the longitude

cosphi = cos(phi);
sinphi = sin(phi);
coslat=cos(lat*pi/180);
sinlat=sin(lat*pi/180);
Bx = Br*coslat*cosphi(1) + Bt*sinlat*cosphi(1) - Bp*sinphi(1);
By = Br*coslat*sinphi(1) + Bt*sinlat*sinphi(1) - Bp*cosphi(1);
Bz = Br*sinlat - Bt*coslat;

B_ECEF = [Bx,By,Bz]';
end

