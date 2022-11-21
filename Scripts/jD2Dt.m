function Date = jD2Dt(jD) 
% Converts time in [yyyy mm dd HH MM secs] to julian date from 1 jan -4713


%%%%%%%%%%%%%%%              WARNING           %%%%%%%%%%%%%%%%%%%%%%%
% This calculation of the Date is verified only for Julian Day after %
%                           2 299 161                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


J = jD + 0.5;

Z = floor(J);  %  Number of day since 1 Jan -4712
F = J - Z;     %  Factional part of the day

if Z < 2299161 % Before 15 oct 1582
    S = Z;
else
    alpha = fix((Z - 1867216.25) / 36524.25);
    S = Z + 1 + alpha - fix(alpha / 4);
end
B = S + 1524;
C = fix((B - 122.1) / 365.25);
D = fix(365.25 * C);
E = fix((B - D) / 30.6001);

% Day of the month
Q = B - D - fix(30.6001 * E); 

% Month of the year
if E < 14
    M = E - 1;
elseif E == 14 || E == 15
    M = E - 13;
else
    error('Error on the calculation of the month of the year');
end

% Year
if M > 2
    A = C - 4716;
elseif M == 1 || M == 2
    A = C - 4715;
else
    error('Error on the calculation of the Year');
end

% Hour
h = round(24 * F);

% Minute
m = round(1440 * (F - h / 24));
if m < 0
    m = m + 60;
    h = h - 1;
end

% Seconde
s = round(86400 * (F - h / 24 - m / 1440));
if s < 0
    s = s + 60;
    m = m - 1;
end

Date = [A M Q h m s];