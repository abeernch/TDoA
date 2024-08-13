function [xout,yout] = coord_tfm(xin,yin,theta,x_offset,y_offset,rdoa)
% Program to rotate and translate x,y values from x",y" to x,y space.
% Written to plot hyperbolas for time of arrival code.
% theta value assumed to be in degrees.
% rotation matrix
xfm = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% Make x,y values into a column vector
r_in = [xin; yin];

%% Rotate
r_out = xfm*r_in;
x = r_out(1,:); y = r_out(2,:);

%% Apply offset to the coordinates. The offset is the half of range between
% the slave and reference sensor.
% if x_offset < 0 && (s1(1) - s2(1))>=0
if rdoa >=0
    xout = x - x_offset;
% elseif x_offset >= 0 && (s1(1)-s2(1)) <= 0
elseif rdoa<0
    xout = x + x_offset;
end
% if y_offset < 0 && (s1(2) - s2(2))>=0 
if rdoa >=0
    yout = y - y_offset;
% elseif y_offset >= 0 && (s1(2) - s2(2))<=0 
elseif rdoa<0
    yout = y + y_offset;
end
end