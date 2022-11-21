function [r, v] = PVOrb( el, t )

%--------------------------------------------------------------------------
%   Generate an orbit by propagating Keplerian elements.
%   Note this function is valid for only ellispoid and circular orbits
%--------------------------------------------------------------------------
%   Form:
%   [r, v, t] = PVOrb( el, t )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   el         Elements vector [a,i,W,w,e,M]
%   t          Times from 0 to ° (sec)

%   -------
%   Outputs
%   -------
%   r          Position vectors in eci (km)for times t
%   v          Velocity vectors in eci(km/sec)for times t

a = el(1);  % semi-major axis
i = el(2);  % Orbit inclination       
L = el(3);  % longitude of the ascending node(rad)
w = el(4);  % argument of perigee(rad)
e = el(5);  % eccentricity

mu = 3.98600436e5;

wo = sqrt(mu/abs(a)^3); % orbit rate

% Transforms from the perifocal frame to the inertial frame
%----------------------------------------------------------
ci = cos(i); cw = cos(w); cL = cos(L);
si = sin(i); sw = sin(w); sL = sin(L);

c = [ cL*cw-sL*sw*ci,-cL*sw-sL*cw*ci, sL*si;...
    sL*cw+cL*sw*ci,-sL*sw+cL*cw*ci,-cL*si;...
    sw*si,          cw*si,    ci];

M      = wo*(t - el(7)) + el(6);

ee = DupVec(e,length(M))';

% only for ellipsoid orbits
k = find(ee < 1);
if( ~isempty(k) )
    eccAnomX(k) = M2EEA(ee(k),M(k));
end

E      = eccAnomX; % eccentric anamoly

if( length(e) == 1 )
    ec = DupVec(e,length(E))';
end

nuX = zeros(size(E));

k = find( ec < 1 );
nuX(k) = 2*atan(sqrt((1+ec(k))./(1-ec(k))).*tan(0.5*E(k)));

theta  = nuX;
cTheta = cos( theta );
sTheta = sin( theta );


aV    = DupVec(a',length(e));
eV    = 1 - e.^2;
i     = find(e == 1);
eV(i) = 2*ones(size(i));
pX    = aV.*DupVec(eV,length(a));
p      = pX;


rMag   = p./(1 + e*cTheta);

r     = c*[rMag.*cTheta;rMag.*sTheta;zeros(size(t))];
v     = sqrt(mu/p)*c*[-sTheta;e+cTheta;zeros(size(t))];
end

function eccAnom = M2EEA( ecc, meanAnom )

tol = 1.e-8;

% Ellipse
if any(ecc >= 1),
    error('The eccentricity must be < 1')
end

e=DupVec(ecc,length(meanAnom));

if any( e < 0 | e == 1 )
    error('The eccentricity must be > 0, and not == 1');
else
    eccAno = zeros(size(meanAnom));
    k    = find( meanAnom ~= 0 );
    e    = e(k);
    m    = meanAnom(k);
    i    = find( m > pi );
    if( ~isempty(i) )
        m(i) = -m(i);
    end
    
    eA   = zeros(size(m));
    
    kL = eval( ['find(', 'e<1',')'] );
    kG = eval( ['find(~(' 'e<1''))'] );
    
    
    if( ~isempty(kL) )     % elliptical case.
        sM      = sin(m(kL));
        eA(kL) = m(kL) + e(kL).*sM./(1 - sin(m(kL)+e(kL)) + sM);
    end
    
    if( ~isempty(kG) )     % hyperbolic case.
        sM     = sinh(m(kG)./(e(kG)-1));
        eA(kG) = m(kG).^2 ./ (e(kG).*(e(kG)-1).*sM - m(kG));
    end;
    
    if( ~isempty(i) )
        eA(i) = -eA(i);
    end
    
    eccAno(k) = eA;
    
end
eccAnomX=eccAno;

% Iterate
delta = tol + 1;
n     = 0;
tau   = tol;

while ( max(abs(delta)) > tau )
    dE    	  = (meanAnom - eccAnomX + ecc.*sin(eccAnomX))./ ...
        (1 - ecc.*cos(eccAnomX));
    eccAnomX    = eccAnomX + dE;
    n           = n + 1;
    delta       = norm(abs(dE),'inf');
    tau         = tol*max(norm(eccAnomX,'inf'),1.0);
end

eccAnom = eccAnomX;
end

function y = DupVec( x, n )

[r,c] = size(x);

if( r > c )
    y = x(:,ones(1,n));
else
    y = x(ones(n,1),:);
end
end