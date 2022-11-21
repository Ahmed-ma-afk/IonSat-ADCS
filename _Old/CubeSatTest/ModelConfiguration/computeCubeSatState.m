function initCond = computeCubeSatState(blk, mode, sim_t0, epoch, varargin)
%computeCubeSatState Compute CubeSat orbit trajectory and attitude

%  Copyright 2019 MathWorks, Inc.

persistent savedInitCond

dAT = 0; dUT1 = 0; pm = [0 0]; dCIP = [0 0]; lod = 0;

% Convert julian dates to date vectors
epochVec = datevec(datetime(epoch,'convertfrom', 'juliandate'));
sim_t0Vec = datevec(datetime(sim_t0,'convertfrom', 'juliandate'));

% Conpute transformation between J2000 and MOD ECI frames
R_ijk2j2000 = dcmIJK2J2000(epochVec, dAT)';

switch mode
    case 'Keplerian Orbital Elements'
        % Parse inputs
        a = varargin{1};
        ecc = varargin{2};
        incl = varargin{3};
        RAAN = varargin{4};
        argp = varargin{5};
        nu = varargin{6};
        truelon = varargin{7};
        arglat = varargin{8};
        lonper = varargin{9};
        euler = varargin{10};
        pqr = varargin{11};
        if nargin > 15
            dAT = varargin{12};
            dUT1 = varargin{13};
            pm = varargin{14};
            dCIP = varargin{15};
            lod = varargin{16};
        end
        
        % Kepler to RV_eci
        small = 1e-12;
        if ( ecc < small )
            if incl < small || abs(incl-pi)< small % circular equatorial orbit
                [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, 0, 0,...
                    0, 'truelon', truelon);
            else % circular inclined
                [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, RAAN, 0,...
                    0, 'arglat', arglat);
            end
        else
            if ( ( incl<small) || (abs(incl-pi)<small) ) % elliptical equatorial
                [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, 0, 0,...
                    nu, 'lonper', lonper);
            else % elliptical inclined
                [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, RAAN, argp,...
                    nu);
            end
        end
        
        % RV_tod to RV_j2000
        r_j2000 = R_ijk2j2000*r_ijk;
        v_j2000 = R_ijk2j2000*v_ijk;
        
        % RV_j2000 to RV_ecef
        [r_ecef, v_ecef] = eci2ecef(sim_t0Vec, r_j2000, v_j2000, 'dAT', dAT,...
            'dUT1', dUT1, 'pm', pm, 'dCIP', dCIP, 'lod', lod);
        
        % R_ecef to R_lla
        lla = ecef2lla(r_ecef(:)');
        
        % V_ecef to V_ned
        v_ned = dcmecef2ned(lla(1), lla(2))*v_ecef(:);
        
        % V_ned to V_body
        uvw = angle2dcm(euler(3)*pi/180, euler(2)*pi/180, euler(1)*pi/180, 'ZYX')*v_ned;
        
    case 'ECI Position and Velocity'
        % Parse inputs
        r_ijk = varargin{1};
        v_ijk = varargin{2};
        euler = varargin{3};
        pqr = varargin{4};
        if nargin > 8
            dAT = varargin{5};
            dUT1 = varargin{6};
            pm = varargin{7};
            dCIP = varargin{8};
            lod = varargin{9};
        end
        
        % RV_eci to Kepler
        [a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] =...
            ijk2keplerian(r_ijk,v_ijk);
        
        % RV_tod to RV_j2000
        r_j2000 = R_ijk2j2000*r_ijk(:);
        v_j2000 = R_ijk2j2000*v_ijk(:);
        
        % RV_j2000 to RV_ecef
        [r_ecef, v_ecef] = eci2ecef(sim_t0Vec, r_j2000, v_j2000, 'dAT', dAT,...
            'dUT1', dUT1, 'pm', pm, 'dCIP', dCIP, 'lod', lod);
        
        % R_ecef to R_lla
        lla = ecef2lla(r_ecef(:)');
        
        % V_ecef to V_ned
        v_ned = dcmecef2ned(lla(1), lla(2))*v_ecef(:);
        
        % V_ned to V_body
        uvw = angle2dcm(euler(3)*pi/180, euler(2)*pi/180, euler(1)*pi/180, 'ZYX')*v_ned;
        
    case 'ECEF Position and Velocity'
        % Parse inputs
        r_ecef = varargin{1};
        v_ecef = varargin{2};
        euler = varargin{3};
        pqr = varargin{4};
        if nargin > 8
            dAT = varargin{5};
            dUT1 = varargin{6};
            pm = varargin{7};
            dCIP = varargin{8};
            lod = varargin{9};
        end
        
        % RV_ecef to RV_j2000
        [r_j2000, v_j2000] = ecef2eci(sim_t0Vec, r_ecef, v_ecef, 'dAT', dAT,...
            'dUT1', dUT1, 'pm', pm, 'dCIP', dCIP, 'lod', lod);
        
        % RV_j2000 to RV_tod
        r_ijk = R_ijk2j2000*r_j2000;
        v_ijk = R_ijk2j2000*v_j2000;
        
        % RV_eci to Kepler
        [a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] =...
            ijk2keplerian(r_ijk,v_ijk);
                
        % R_ecef to R_lla
        lla = ecef2lla(r_ecef(:)');
        
        % V_ecef to V_ned
        v_ned = dcmecef2ned(lla(1), lla(2))*v_ecef(:);
        
        % V_ned to V_body
        uvw = angle2dcm(euler(3)*pi/180, euler(2)*pi/180, euler(1)*pi/180, 'ZYX')*v_ned;
        
    case 'Geodetic LatLonAlt and Velocity in NED'
        % Parse inputs
        lla = varargin{1};
        v_ned = varargin{2};
        euler = varargin{3};
        pqr = varargin{4};
        if nargin > 8
            dAT = varargin{5};
            dUT1 = varargin{6};
            pm = varargin{7};
            dCIP = varargin{8};
            lod = varargin{9};
        end
        
        % V_lla to V_ecef
        r_ecef = lla2ecef(lla(:)');
        
        % V_ned to V_body 
        uvw = angle2dcm(euler(3)*pi/180, euler(2)*pi/180, euler(1)*pi/180, 'ZYX')*v_ned(:);
        
        % V_ned to V_ecef
        v_ecef = dcmecef2ned(lla(1), lla(2))'*v_ned(:);
        
        % RV_ecef to RV_j2000
        [r_j2000, v_j2000] = ecef2eci(sim_t0Vec, r_ecef, v_ecef, 'dAT', dAT,...
            'dUT1', dUT1, 'pm', pm, 'dCIP', dCIP, 'lod', lod);
        
        % RV_j2000 to RV_tod
        r_ijk = R_ijk2j2000*r_j2000;
        v_ijk = R_ijk2j2000*v_j2000;
        
        % RV_eci to Kepler
        [a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] =...
            ijk2keplerian(r_ijk,v_ijk);
end

[~, thGAST] = greenwichSRT(sim_t0, dUT1, dAT);

% Update initial conditions in base workspace
initCond.simStartDate.JD = sim_t0;
initCond.simStartDate.dateVector = sim_t0Vec;
initCond.CoordEpoch.JD = epoch;
initCond.CoordEpoch.dateVector = epochVec;
initCond.OrbitalElements.semiMajorAxis = a;
initCond.OrbitalElements.eccentricity = ecc;
initCond.OrbitalElements.inclination = incl;
initCond.OrbitalElements.RAAN = RAAN;
initCond.OrbitalElements.argumentOfPerigee = argp;
initCond.OrbitalElements.trueAnomoly = nu;
initCond.OrbitalElements.trueLongitude = truelon;
initCond.OrbitalElements.argumentOfLatititude = arglat;
initCond.OrbitalElements.longitudeOfPerigee = lonper;
initCond.r_eci = r_ijk(:)';
initCond.v_eci = v_ijk(:)';
initCond.r_ecef = r_ecef(:)';
initCond.v_ecef = v_ecef(:)';
initCond.lla = lla(:)';
initCond.v_ned = v_ned(:)';
initCond.uvw = uvw(:)';
initCond.euler = euler(:)';
initCond.pqr = pqr(:)';
initCond.EarthProps.dAT = dAT;
initCond.EarthProps.dUT1 = dUT1;
initCond.EarthProps.pm = pm;
initCond.EarthProps.dCIP = dCIP;
initCond.EarthProps.lod = lod;
initCond.EarthProps.LG = thGAST;
assignin('base','initCond',initCond);

% Sync initial conditions block with base workspace if values have changed
if isempty(savedInitCond) || ~isequal(initCond,savedInitCond)
    syncICwithWS(blk, initCond);
    savedInitCond = initCond;
end

end

function dcm = dcmIJK2J2000(epochVec, dAT)
%Compute tranformation from ECI with mean equinox at epoch to J2000

% Seconds for UTC
ssTT = epochVec(end) + dAT + 32.184;
% Julian date for terrestrial time
jdTT = mjuliandate(epochVec(1),epochVec(2),epochVec(3),epochVec(4),epochVec(5),ssTT);
% Number of Julian centuries since J2000 for terrestrial time.
tTT = (jdTT - 51544.5)/36525;
tTT2 = tTT.*tTT;
tTT3 = tTT2.*tTT;
% Zeta, theta and z represent the combined effects of general precession
zeta = convang((2306.2181*tTT + 0.30188*tTT2 + 0.017998*tTT3)/3600,'deg','rad');
theta = convang((2004.3109*tTT - 0.42665*tTT2 - 0.041833*tTT3)/3600,'deg','rad');
z = convang((2306.2181*tTT + 1.09468*tTT2 + 0.018203*tTT3)/3600,'deg','rad');
% ECI vector with mean equinox at epoch to ECI vector in J2000
dcm = angle2dcm(-zeta,theta,-z,'ZYZ')';
end

function syncICwithWS(blk, initCond)
%syncICwithWS Update modified block parameters in Initial Conditions block

set_param(blk, 'sim_t0' ,num2str(initCond.simStartDate.JD));
set_param(blk, 'epoch' ,num2str(initCond.CoordEpoch.JD));
set_param(blk, 'a', num2str(initCond.OrbitalElements.semiMajorAxis));
set_param(blk, 'ecc', num2str(initCond.OrbitalElements.eccentricity));
set_param(blk, 'incl', num2str(initCond.OrbitalElements.inclination));
set_param(blk, 'omega', num2str(initCond.OrbitalElements.RAAN));
set_param(blk, 'argp', num2str(initCond.OrbitalElements.argumentOfPerigee));
set_param(blk, 'nu', num2str(initCond.OrbitalElements.trueAnomoly));
set_param(blk, 'truelon', num2str(initCond.OrbitalElements.trueLongitude));
set_param(blk, 'arglat', num2str(initCond.OrbitalElements.argumentOfLatititude));
set_param(blk, 'lonper', num2str(initCond.OrbitalElements.longitudeOfPerigee));
set_param(blk, 'r_eci', regexprep(['[' num2str(initCond.r_eci(:)') ']'],' +',' '));
set_param(blk, 'v_eci', regexprep(['[' num2str(initCond.v_eci(:)') ']'],' +',' '));
set_param(blk, 'r_ecef', regexprep(['[' num2str(initCond.r_ecef(:)') ']'],' +',' '));
set_param(blk, 'v_ecef', regexprep(['[' num2str(initCond.v_ecef(:)') ']'],' +',' '));
set_param(blk, 'lla', regexprep(['[' num2str(initCond.lla(:)') ']'],' +',' '));
set_param(blk, 'v_ned', regexprep(['[' num2str(initCond.v_ned(:)') ']'],' +',' '));
set_param(blk, 'euler', regexprep(['[' num2str(initCond.euler(:)') ']'],' +',' '));
set_param(blk, 'pqr', regexprep(['[' num2str(initCond.pqr(:)') ']'],' +',' '));
set_param(blk, 'dAT', num2str(initCond.EarthProps.dAT));
set_param(blk, 'dUT1', num2str(initCond.EarthProps.dUT1));
set_param(blk, 'pm', regexprep(['[' num2str(initCond.EarthProps.pm(:)') ']'],' +',' '));
set_param(blk, 'dCIP', regexprep(['[' num2str(initCond.EarthProps.dCIP(:)') ']'],' +',' '));
set_param(blk, 'lod', num2str(initCond.EarthProps.lod));
end

