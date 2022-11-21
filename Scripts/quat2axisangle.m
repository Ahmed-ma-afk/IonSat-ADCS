function [y] = quat2axisangle(u)
%UNTITLED Summary of this function goes here
%   Inputs: quaternion
%   Outputs: axis angle [x y z theta (rad)]
   if (u(1) > 1) 
       u=u/norm(u); %if w>1 acos and sqrt will produce errors, this cant happen if quaternion is normalised
   end
qw=u(1);
qx=u(2);
qy=u(3);
qz=u(4);   
   
   angle = 2 * acos(qw);
   if angle > pi
       angle=angle-2*pi;
   end
   s = sqrt(1-qw*qw); % assuming quaternion normalised then w is less than 1, so term always positive.
   if (s < 0.001) 
       % test to avoid divide by zero, s is always positive due to sqrt
        % if s close to zero then direction of axis not important
        x = qx; % if it is important that axis is normalised then replace with x=1; y=z=0;
        y = qy;
        z = qz;
   else
       x = qx / s; % normalise axis
       y = qy / s;
       z = qz / s;
   end

y = [x;y;z;angle];
end

