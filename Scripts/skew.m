function [skew_matrix] = skew(u)
%SKEW Summary of this function goes here
%   Detailed explanation goes here
u1=u(1);
u2=u(2);
u3=u(3);
skew_matrix = [  0 -u3  u2 ;...
                u3   0 -u1 ;...
               -u2  u1   0 ];
end

