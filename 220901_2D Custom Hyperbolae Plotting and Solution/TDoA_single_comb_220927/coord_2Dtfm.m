function [xout,yout] = coord_2Dtfm(xin,yin,theta,x_offset,y_offset,rdoa)
%% Performs 2-D coordinate transformation on a hyperbola
% This is a nested function that is called within the isochrone_func. 
% - Once a generic hyperbola is created, the necassary coordinate transformation from body to fixed frame
% reference is completed here  w.r.t the reference sensor. 
% - The transformation angle (theta), adjusted offset center of the hyperbola and RDOA are input. 
% - The Rotation matrix is calculated and matrix multiplied with the coordinates of the hyperbolae
% - Based on the sign of RDOA, the calculated offset is either added (-ive rdoa) or subtracted (+ive rdoa) from the rotated points 

%% Compute the roatation matrix
xfm = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

%% Make x,y values into a column vector
r_in = [xin; yin];

%% Rotate
r_out = xfm*r_in;
x = r_out(1,:); y = r_out(2,:);

%% Apply offset to the coordinates. 
% The offset is the half of range between the slave and reference sensor.
% if x_offset < 0 && (s1(1) - s2(1))>=0
if rdoa >= 0
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

%% Change log:
% 1. Created: 220926