%% Plotting 2 sheet hyperboloid
% Essentially a rotation of the 2D hyperbola along the axis that bisects it

% Script testing creation of the quadric surface (elliptical 2-sheet
% hyperboloid) for 3D TDoA algorithm.

% Test 3D coordinate transformation

clear;clc

%% Basic scenario with sensor and source pos
sensor_pos = [0e3 5e3;
              0e3 2.5e3;
              0e3 1e3];

source_pos = [-2e3;
               0.5e3;
               1.5e3];

ranges = vecnorm(source_pos - sensor_pos);
rdoa = ranges(1) - ranges(2);

%% Hyperboloid configuration parameters 
% Important for coordinate transformation
% a = 1; % 
% b = 0.1;
% c = rdoa; %???????????????????????????

% a = abs(rdoa)/2;
% d = norm(sensor_pos(:,1) - sensor_pos(:,2))/2;
% b = sqrt(d^2 - a^2);
% c = norm(sensor_pos(3,2) - sensor_pos(3,1)) %sqrt(a^2 + b^2);
% Toggle between hyperboloid or its projection in 2D

%% Point symmetry about z-axis
a = 0.1; b = 0.1; c = 1%abs(rdoa)/2e3;

D2 = 0;
[min_s, max_s] = deal(-20e3,20e3);
if D2 == 0
    [X,Y] = meshgrid(min_s:0.2*((max_s-min_s)/10):max_s); % 51 point grid
else
    [X] = meshgrid(min_s:0.2*((max_s-min_s)/10):max_s);
    Y = zeros(size(X));
end
% x = linspace(-20e3,20e3,20e3);
% y = linspace(-20e3,20e3,20e3);
Z = (c*sqrt(a^2*b^2 + a^2*Y.^2 + b^2*X.^2))/(a*b);
% Z = c*sqrt((x/a).^2 + (y/b).^2 + 1);
% figure
% surf(X,Y,-Z)
% hold on
% surf(X,Y,Z)
% hold off
% grid on
% xlim([min_s max_s]);ylim([min_s max_s]);zlim([min_s max_s])
% % shading interp
% axis  vis3d 
%% Point symmetry about x-axis
% if D2 == 0
%     [X,Z] = meshgrid(min:0.2*((max-min)/10):max); % 51 point grid
% else
%     [X] = meshgrid(min:0.2*((max-min)/10):max); % 51 point grid
%     Z = zeros(size(X));
% end
% 
% Y = (b/(a*c))*sqrt(a^2*Z.^2 + c^2*X.^2 + a^2*c^2);
% 
% figure
% surf(X,-Y,Z)
% hold on
% surf(X,Y,Z)
% hold off
% grid on
% % shading interp
% xlim([min max]);ylim([min max]);zlim([min max])
% % axis  vis3d 

%% Test 3D coordinate transformation
x = X(:); y = Y(:); z = Z(:);
center = (sensor_pos(:,2) - sensor_pos(:,1))/2;
% Azimuth and elevation angles with respecto the fixed frame (reference sensor)
az = atan2d((sensor_pos(2,2) - sensor_pos(2,1)),((sensor_pos(1,2) - sensor_pos(1,1)))) + 0;
el = atan2d((sensor_pos(3,2) - sensor_pos(3,1)),(sqrt((sensor_pos(1,2) - sensor_pos(1,1))^2 + (sensor_pos(2,2) - sensor_pos(2,1))^2)))+ 90;

for i = 1:length(x)
    [xout,yout,zout] = coord3D_tfm(x(i),y(i),z(i),az,el,center(1),center(2),center(3),rdoa);
    x_trf(i) = xout; 
    y_trf(i) = yout;
    z_trf(i) = zout;
end

if rdoa < 0
    x_trf = x_trf + sensor_pos(1,1);
    y_trf = y_trf + sensor_pos(2,1);
    z_trf = z_trf + sensor_pos(3,1);
elseif rdoa >= 0
    x_trf = -x_trf + sensor_pos(1,1);
    y_trf = -y_trf + sensor_pos(2,1);
    z_trf = -z_trf + sensor_pos(3,1);
end

x_trf = reshape(x_trf,size(X,1),size(X,2));
y_trf = reshape(y_trf,size(Y,1),size(Y,2));
z_trf = reshape(z_trf,size(Z,1),size(Z,2));

figure; s =surf(x_trf,y_trf,z_trf); s.FaceAlpha = 0.25;
hold on
surf(-x_trf ,-y_trf,-z_trf,z_trf,"FaceAlpha",0.25)
hold on
scatter3(sensor_pos(1,:),sensor_pos(2,:),sensor_pos(3,:),50,"blue","filled","o")
scatter3(source_pos(1,:),source_pos(2,:),source_pos(3,:),50,"red","filled","^")
axis([min_s max_s min_s max_s min_s max_s]) % specify axis limits


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