 clear;clc
syms x1 x2 x3 x y1 y2 y3 y t1 t2 t3 t c
% eqns = [c^2*(t1^2 - 2*t1*t - t2^2 +2*t2*t)- x1^2 + x2^2 - y1^2 - y2^2 + 2*x1*x - 2*x2*x + 2*y1*y - 2*y2*y==0,],

eqn1 = sqrt((x1-x)^2 + (y1-y)^2) - sqrt((x2-x)^2 + (y2-y)^2)  == c*(t1-t2);
eqn2 = sqrt((x1-x)^2 + (y1-y)^2) - sqrt((x3-x)^2 + (y3-y)^2)  == c*(t1-t3);

% sol = solve([eqn1 eqn2],[x y]);



%% Define sensor and source position
[x1,y1]       = deal(0,0);
[x2,y2]       = deal(5e3,5e3);
[x3,y3]       = deal(5e3,0);
[x_act,y_act] = deal(25e3, 1e3);
c = 3e8;

%% Calculate sensor-source range
ranges = rangeangle([x1 x2 x3;y1 y2 y3; 0 0 0],[x_act;y_act;0])

%% Convert ranges to TOA
toa = range2time(ranges,c)*0.5;
t1 = toa(1);
t2 = toa(2);
t3 = toa(3);


%% Evaluate equations
s1 = eval(eqn1);
s2 = eval(eqn2);
sol = solve([s1 s2],[x y]);
%% Extract soln
xsol = sol.x;
ysol = sol.y;
x_est = (eval(xsol))
y_est = (eval(ysol))