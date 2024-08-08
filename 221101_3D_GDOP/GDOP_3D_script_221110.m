%% Script for 3D GDOP computation

% 221110
%% Initialize inputs to the function

% Sensor Position
sensor_pos = [00e3 20e3 -20e3 00e3;
              00e3 15e3 15e3 -15e3;
              00e3 0.0e3 0.25e3 0.1e3];

% sensor_pos = [0	3116.05483775860	-11496.8490106353;
%               0	-7246.75408749039	-21.4179927785758;
%               0	-75.8919118246145	46.6497395972486];

% Sensor timing error
timing_err = 50e-9;

% Grid resolution (in meters) along x, y and z axis
grd_res_xyz = [1000; 1000; 100];

% Axis limits
x_span = [-400e3 400e3];
y_span = [-400e3 400e3];
z_span = [-10e3 10e3];

res = [(sum(abs(x_span)))/grd_res_xyz(1);(sum(abs(y_span)))/grd_res_xyz(2);(sum(abs(z_span)))/grd_res_xyz(3)];

%% Function
[gdop_xy,gdop_xz,gdop_yz] = GDOP_3dfn(sensor_pos,timing_err,res,x_span,y_span,z_span);