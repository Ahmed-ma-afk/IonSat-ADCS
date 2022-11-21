function jD = Dt2jD(Date) 
% Converts time in [yyyy mm dd HH MM secs] to julian date from 1 jan -4713


%%%%%%%%%%%%%%%              WARNING           %%%%%%%%%%%%%%%%%%%%%%%
% This calculation of the Julian Day is verified only for date after %
%                          15 Otc 1582                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = Date(1);    % Year
M = Date(2);    % Month
D = Date(3);    % Day
% Fractional part of the Day
F = (Date(4) - 12 + (Date(5) + Date(6) / 60) / 60) / 24;

if( M <= 2 )
  M = M + 12;
  A = A - 1;
end
S=fix(A / 100);
B=2 - S + fix(S / 4);
jD = fix(365.25 * (A + 4716)) + B + fix(30.6001 * (M + 1)) + D + F -1524;