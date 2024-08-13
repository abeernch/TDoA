
%% Coord trm function
function [xout,yout,zout] = coord3D_tfm(xin,yin,zin,az,el,x_offset,y_offset,z_offset,rdoa)
% Program to rotate and translate x,y,z values from x",y",z" to x,y,z space.
% Written to plot hyperbooloids for 3D time of arrival code.
% theta value assumed to be in degrees.

%% 3D rotation matrix
xfm = [cosd(az)*cosd(el) -sind(az) -cosd(az)*sin(el);
       sind(az)*cosd(el) cosd(az)  -sind(az)*sind(el);
       sind(el)             0            cosd(el)];

%% Make x,y,z values into a column vector
r_in = [xin; yin; zin];

%% Rotate
r_out = xfm*r_in;
x = r_out(1,:); y = r_out(2,:); z = r_out(3,:);

%% Apply offset to the coordinates. The offset is the half of range between
% the slave and reference sensor.

if rdoa >= 0
    xout = x - x_offset;
    yout = y -  y_offset;
    zout = z -  z_offset;
elseif rdoa < 0
    xout = x + x_offset;
    yout = y +  y_offset;
    zout = z +  z_offset;
end
end