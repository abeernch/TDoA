function iso = draw_2Disochrone(s1,s2,rdoa,x_lim,points)
%% CREATES GENERIC HYPERBOLA AND COMPUTES PARAMETERS FOR 2D COORD TFM
% -This function creates hyperbolae in the opening up/down configuration. 
% -The inputs to the function are:
    % 1. Sensor network geometry parameters, 
    % 2. 4-quad angle of the slave sensor w.r.t the reference sensor, 
    % 3. the max orthogonal distance of hyperbola, 
    % 4. no of points to plot and 
    % 5. RDOA with respect to the reference sensor. 

% -This function is called iteratively for n tdoa pairs to plot n hyperbolae.

% Author: Abeer Nasir Chaudhry®
%% Compute components of the hyperbolic geometry

a = abs(rdoa)/2;          % Calculate distance from vertix to center of hyperbola
d = norm(s1-s2)/2;        % Calculate distance from foci to center of hyperbola (half the distance between two sensors (foci))
b = sqrt(d^2 - a^2);      % Conjugate axis parameter to determine the hyperbola
center = ((s2-s1)/2);     % Find the center offset points for a generic hyperbola w.r.t tghe sensors' position

%% Calculate angle between slave and reference sensor 
theta = atan2d((s2(2) - s1(2)),((s2(1) - s1(1)))) + 90;
%% Create generic Hyperbola 
xmax = x_lim;                         % Max orthogonal distance of the hyperbola   
x = linspace(-xmax,xmax,points);      % array of x values for plot based on 'points'
y = (a/b).*(sqrt(x.^2 + b^2));        % corresponding y values

%% Apply coord transformation
for i = 1:length(x)
    [xout,yout] = coord_2Dtfm(x(i),y(i),theta,center(1),center(2),rdoa);
    x_trf(i) = xout; 
    y_trf(i) = yout;
end

%% Choose hyperbola leaf
% And perform final offset
if rdoa < 0
    iso = [x_trf;y_trf] + s1(1:2);
elseif rdoa >= 0
    iso = [-x_trf;-y_trf] + s1(1:2);
end
end

% LOG:
% 1. Date created (updated): 220926 (221002)