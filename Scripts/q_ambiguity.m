function [y] = q_ambiguity(u)
%Q_AMBIGUITY Summary of this function goes here
%   Detailed explanation goes here
q0=u(1);
if abs(q0)<0.5
    y=-u;
else
    y=u;
end
    
end

