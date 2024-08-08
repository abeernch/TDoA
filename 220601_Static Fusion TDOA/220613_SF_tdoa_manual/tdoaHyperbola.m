function varargout = tdoaHyperbola(tdoaDetection,type)

emissionSpeed = physconst('Lightspeed');
timeScale = 1e9;

% Location of receiver in global coordinates
p1 = tdoaDetection.MeasurementParameters(1);
r1Pos = p1.OriginPosition;

% Location of receiver 2 in global coordinates
p2 = tdoaDetection.MeasurementParameters(2);
r2Pos = p2.OriginPosition;

theta = linspace(-pi/2,pi/2,100);
if strcmpi(type,'2D')
    phi = 0;
else
    phi = linspace(-pi,pi,50);
end

[THETA, PHI] = meshgrid(theta,phi);

% Compute hyperbola quantities
D = norm(r2Pos - r1Pos)/2;
c = tdoaDetection.Measurement*emissionSpeed/timeScale;

% Hyperbola equation in a frame where sensors lie on x axis and are located
% at -D and D on x axis
xC = -c./cos(THETA)./2;
yC = sqrt(4*D^2 - c^2).*tan(THETA).*cos(PHI)./2;
zC = sqrt(4*D^2 - c^2).*tan(THETA).*sin(PHI)./2;

% Translate and rotate to the scenario frame
r0 = (r1Pos + r2Pos)/2;
a = [1;0;0];
b = (r1Pos - r2Pos);
b = b/norm(b);
v = cross(a,b);
s = norm(v);
c = dot(a,b);
V = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
if abs(s) > 0
    R = eye(3) + V + V^2*(1 - c)/s^2;
else
    R = eye(3);
end

x = R(1,1).*xC + R(1,2).*yC + R(1,3).*zC;
y = R(2,1).*xC + R(2,2).*yC + R(2,3).*zC;
z = R(3,1).*xC + R(3,2).*yC + R(3,3).*zC;
x = x + r0(1);
y = y + r0(2);
z = z + r0(3);

if strcmpi(type,'2D')
    x = x(:);
    y = y(:);
    z = z(:);
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
else
    varargout{1} = surf2patch(x,y,z);
end

end