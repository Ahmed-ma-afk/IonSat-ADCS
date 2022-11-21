function [eclipse] = eclipse_check(u)
%Checks if the satellite is in eclipse
% Inputs:
% u(1:3)= position of the satellite in ECI
% u(4:6)= direction of the Sun in ECI
%   
% Output: out=1: S/C in eclipse 
%         out=0: S/C in illumination

%Defining the frames xyz=inertial XYZ=z pointing to sun
Z=u(4:6);
if Z==[0;0;0]
    X=[1;0;0];
else 
    X=[-Z(2);Z(1);0];
end

X=X/norm(X);
Y=cross(Z,X);

R_i2s=[ X'*[1;0;0] , Y'*[1;0;0] , Z'*[1;0;0] ;...
        X'*[0;1;0] , Y'*[0;1;0] , Z'*[0;1;0] ;...
        X'*[0;0;1] , Y'*[0;0;1] , Z'*[0;0;1] ];
%position of the satellite in new frame XYZ=z pointing to sun    
r_s=R_i2s'*u(1:3);
a=r_s(1);
b=r_s(2);
c=r_s(3);
%check if satellite is in shadow: shadow=0
if c<0
    if ((b^2)+(a^2))<4.0678884e+13
        eclipse=0;
    else
        eclipse=1;
    end
else
    eclipse=1;
end