function [xunit, yunit] = circle(x,y,r,pts)
th = linspace(-pi/2,pi/2,pts);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
