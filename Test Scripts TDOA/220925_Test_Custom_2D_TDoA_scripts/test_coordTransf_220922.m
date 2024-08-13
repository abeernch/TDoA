% Program to rotate and translate x,y values from x",y" to x,y space.
% Written to plot hyperbolas for time of arrival code.
% theta value assumed to be in radians.
% rotation matrix
xfm 
% make x,y values into a column vector
r_in = [xin; yin];
% rotate
r_out = xfm*r_in;
x = r_out(1,:); y = r_out(2,:);
xout = x;
yout = y;