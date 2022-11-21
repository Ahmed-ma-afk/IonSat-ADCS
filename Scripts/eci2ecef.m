function [DCM_eci2ecef] = eci2ecef(UTC)
%Gets the Direction Cosine Matrix for conversion between ECI and ECEF
%   Input: UTC time (6x1)
%   Output: DCM matrix (3x3)
year =  UTC(1);
mes =   UTC(2);
dia   = UTC(3);
hora  = UTC(4);
minuto= UTC(5);
segundo=UTC(6);

%1: get modified julian day (mjd)
mjd = 367*year + dia - 712269 + fix(275*mes/9)- fix(7*(year+fix((mes+9)/12))/4);
%2: get fraction of the day (dfra)
dfra = segundo + 60*(minuto + 60*hora);
%3: get Greenwich Sidereal Time (GWST)
tsj = (mjd-18262.5)/36525;
tsgo = (24110.54841 + (8640184.812866 + 9.3104e-2*tsj - 6.2e-6*tsj*tsj)*tsj)*pi/43200;
tetp = 7.292116e-5;		% Earth angular velocity (rad/s)
gwst = mod(tsgo + dfra*tetp, 2*pi);
%4: get DCM from ECI to ECEF
coan    =  cos(gwst);
sian    =  sin(gwst);

DCM_eci2ecef =  [coan, sian, 0; -sian, coan, 0; 0, 0, 1];
end

