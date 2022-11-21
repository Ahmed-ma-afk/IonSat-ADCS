function [gh] = loadigrfcoefsim(time)

% LOADIGRFCOEFS Load coefficients used in IGRF model.
% 
% Usage: [G, H] = LOADIGRFCOEFS(TIME) or GH = LOADIGRFCOEFS(TIME)
% 
% Loads the coefficients used in the IGRF model at time TIME in MATLAB
% serial date number format and performs the necessary interpolation. If
% two output arguments are requested, this returns the properly
% interpolated matrices G and H from igrfcoefs.mat. If just one output is
% requested, the proper coefficient vector GH from igrfcoefs.mat is
% returned.
% 
% If this function cannot find a file called igrfcoefs.mat in the MATLAB
% path, it will try to create it by calling GETIGRFCOEFS.
% 
% Inputs:
%   -TIME: Time to load coefficients either in MATLAB serial date number
%   format or a string that can be converted into MATLAB serial date number
%   format using DATENUM with no format specified (see documentation of
%   DATENUM for more information).
% 
% Outputs:
%   -G: g coefficients matrix (with n going down the rows, m along the
%   columns) interpolated as necessary for the input TIME.
%   -H: h coefficients matrix (with n going down the rows, m along the
%   columns) interpolated as necessary for the input TIME.
%   -GH: g and h coefficient vector formatted as:
%   [g(n=1,m=0) g(n=1,m=1) h(n=1,m=1) g(n=2,m=0) g(n=2,m=1) h(n=2,m=1) ...]
% 
% See also: IGRF, GETIGRFCOEFS.
% Convert time to a datenumber if it is a string.

% Make sure time has only one element.
if numel(time) > 1
    error('loadigrfcoefs:timeInputInvalid', ['The input TIME can only ' ...
        'have one element']);
end
% Convert time to fractional years.
timevec = datevec(time);
annee=datenum([timevec(1) 1 1]);
if (~mod(timevec(1),4) && mod(timevec(1),100)) || (~mod(timevec(1),400))
    time=timevec(1) + (time - annee)./366;
else
    time=timevec(1) + (time - annee)./365;
end

% Load coefs and years variables.

load IGRFsimCoefs.mat;


% Check validity on time.
years = cell2mat({coefsim.year});
if time < years(1) || time > years(end)
    error('igrf:timeOutOfRange', ['This IGRF is only valid between ' ...
        num2str(years(1)) ' and ' num2str(years(end))]);
end

% Get the nearest epoch that the current time is between.
lastepoch = find(years - time < 0, 1, 'last');
nextepoch = lastepoch + 1;


% Get the coefficients based on the epoch.
lastgh = coefsim(lastepoch).gh;
nextgh = coefsim(nextepoch).gh;

% If one of the coefficient vectors is smaller than the other, enlarge
% the smaller one with 0's.
if length(lastgh) > length(nextgh)
    smalln = length(nextgh);
    nextgh = zeros(size(lastgh));
    nextgh(1:smalln) = coefsim(nextepoch).gh;
elseif length(lastgh) < length(nextgh)
    smalln = length(lastgh);
    lastgh = zeros(size(nextgh));
    lastgh(1:smalln) = coefsim(lastepoch).gh;
end

% Calculate gh using a linear interpolation between the last and next
% epoch.
if coefsim(nextepoch).slope
    ghslope = nextgh;
else
    ghslope = (nextgh - lastgh)/diff(years([lastepoch nextepoch]));
end
gh = lastgh + ghslope*(time - years(lastepoch));
    