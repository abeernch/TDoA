function toa = MeasureTDOA(pos, params)
% This function calculates the expected TDOA measurement given object
% position and measurement parameters of a TDOA detection
timeScale = 1e9;
emissionSpeed = physconst('Lightspeed');
r1 = norm(pos(1:3) - params(1).OriginPosition);
r2 = norm(pos(1:3) - params(2).OriginPosition);
toa = (r1 - r2)/emissionSpeed*timeScale;
end