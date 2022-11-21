function q = Qorb( r, v )

%-------------------------------------------------------------------------------
%   Generate the quaternions that transform from ECI to LVLH (orbital)coordinates.
%   For LVLH coordinates;
%   z is in the -r direction (correction --> x is in the -r direction)
%   y is in the - rxv direction
%   x completes the set(in the direction of velocity)
%-------------------------------------------------------------------------------
%   Form:
%   q = Qorb( r, v )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   r          (3,1) Position vectors in km
%   v          (3,1) Velocity vectors in km
%
%   -------
%   Outputs
%   -------
%   q          (4,1) Quaternions C1 %

y       =  cross( r, v );
yn      =  y/sqrt(sum(y.^2));
x       =  -r;
xn      =  x/sqrt(sum(x.^2));
z       =  cross( x, y );
zn      =  z/sqrt(sum(z.^2));
m       = [xn yn zn];
% q        = zeros(4,1);
% 
% q(1)=0.5*sqrt(1 + m(1,1) + m(2,2) + m(3,3));
% q(2)=sign(m(3,2) - m(2,3))*0.5*sqrt(1 + m(1,1) - m(2,2) - m(3,3));
% q(3)=sign(m(1,3) - m(3,1))*0.5*sqrt(1 - m(1,1) + m(2,2) - m(3,3));
% q(4)=sign(m(2,1) - m(1,2))*0.5*sqrt(1 - m(1,1) - m(2,2) + m(3,3));

q = rot2Quat( m );


function q = rot2Quat(m)
q        = zeros(4,1);

v        = [ trace(m), m(1,1), m(2,2), m(3,3) ];

[~,i] = max(v);

if ( i == 1 ),
    q(1) = sqrt( 1 + v(1) );
    q(2) = (m(3,2) - m(2,3)) / q(1);
    q(3) = (m(1,3) - m(3,1)) / q(1);
    q(4) = (m(2,1) - m(1,2)) / q(1);
    
elseif ( i == 2 ),
    q(2) = sqrt( 1 + 2*m(1,1) - v(1) );
    q(1) = (m(3,2) - m(2,3)) / q(2);
    q(3) = (m(1,2) + m(2,1)) / q(2);
    q(4) = (m(1,3) + m(3,1)) / q(2);
    
elseif ( i == 3 ),
    q(3) = sqrt( 1 + 2*m(2,2) - v(1) );
    q(1) = (m(1,3) - m(3,1)) / q(3);
    q(2) = (m(1,2) + m(2,1)) / q(3);
    q(4) = (m(2,3) + m(3,2)) / q(3);
    
elseif ( i == 4 ),
    q(4) = sqrt( 1 + 2*m(3,3) - v(1) );
    q(1) = (m(2,1) - m(1,2)) / q(4);
    q(2) = (m(1,3) + m(3,1)) / q(4);
    q(3) = (m(2,3) + m(3,2)) / q(4);
end

% Halve to get q and make sure that q(1) is positive
%---------------------------------------------------
if q(1) < 0,
    q =  - 0.5 * q;
else
    q =   0.5 * q;
end
