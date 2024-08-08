function [gdop_xy,gdop_xz,gdop_yz,xy_ax] = GDOP_3dfn(sensor_pos,timingerr,grd_res,x_span,y_span,z_span)

% The function takes the following inputs:
% 1. Sensor Positions (in the following format [x1 x2 x3...xn; 
%                                               y1 y2 y3...yn;
%                                               z1 z2 z3...zn])
% 2. Limits for x and y
% 3. TDOA Accuracy (timing error)
% And computes the Geometric dilution of precision (GDOP) by initializing
% possible source positions over a volume and computing the ratio of error
% in position to the error in range for ach pair of sensors.

% DOP represents the positional precision due to the angular separation
% between the sensors that are used to estimate a source's position.

% Log at end
% Author: Abeer Nasir Chaudhry®

%% Check input arguments and assign default values
if nargin == 0
%     sensor_pos = [00e3 5e3 10e3 15e3;00e3 5e3 10e3 20e3; 1e3 0.50e3 0.250e3 00e3];
    sensor_pos = [00e3 -25e3 25e3 00e3;00e3 25e3 25e3 -25e3; 1e3 0.50e3 0.250e3 00e3];

    timingerr = 30e-9;
    grd_res = 1e3;
    x_span = [-200e3 200e3];
    y_span = [-200e3 200e3];
    z_span = [-200e3 200e3];
elseif nargin < 2
    timingerr = 30e-9;
    grd_res = 1e3;
    x_span = [-200e3 200e3];
    y_span = [-200e3 200e3];
    z_span = [-200e3 200e3];
elseif nargin < 3
    grd_res = 1e3;
    x_span = [-200e3 200e3];
    y_span = [-200e3 200e3];
    z_span = [-200e3 200e3];
elseif nargin < 4
    x_span = [-200e3 200e3];
    y_span = [-200e3 200e3];
    z_span = [-200e3 200e3];
elseif nargin < 5
    y_span = [-200e3 200e3];
    z_span = [-200e3 200e3];
elseif nargin < 6
    z_span = [-200e3 200e3];
end

nSensors = size(sensor_pos,2);

nDim = size(sensor_pos,1);
if nDim ~= 3
    error(['Enter sensor positions in the following format [x1 x2 x3...xn; ' ...
                                                          ' y1 y2 y3...yn; ' ...
                                                          ' z1 z2 z3...zn])'])
end

% Extract plane-wise sensor positions in xy, xz and yz 
pos_xy = sensor_pos(1:2,:);
pos_xz = sensor_pos(1:2:3,:);
pos_yz = sensor_pos(2:3,:);

% Assign custom grid resolutions (no. of points calculated from resolution defined earlier)
res_x = grd_res(1);
res_y = grd_res(2);
res_z = grd_res(3);

%% Initialize functions for plane-wise calculation of GDOP
% XY plane
plane = 'xy';
[GDOPxy,x,y] = gen_GDOPfn(pos_xy,plane,timingerr,res_x,res_y,x_span(1),x_span(2),y_span(1),y_span(2));

% XZ plane
plane = 'xz';
[GDOPxz,x,z] = gen_GDOPfn(pos_xz,plane,timingerr,res_x,res_z,x_span(1),x_span(2),z_span(1),z_span(2));

% XY plane
plane = 'yz';
[GDOPyz,y,z] = gen_GDOPfn(pos_yz,plane,timingerr,res_y,res_z,y_span(1),y_span(2),z_span(1),z_span(2));


%% Plotting results

contourLevelsxy = [.1,0.5,1,2,3,5,10,15,20,30,40];
contourLevels = [0.01:0.01:0.5, 0.5:0.1:2, 4:2:10];

%% Draw Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    XY plane      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
xy_ax=axes;
plot(xy_ax,sensor_pos(1,:)/1e3,sensor_pos(2,:)/1e3,'o','DisplayName','Sensors','LineWidth',1, ...
    'MarkerFaceColor','blue','MarkerEdgeColor','blue');
for j = 1:nSensors
    lbl_sen = sprintf('S_{%1.0d}',j);
    text(xy_ax,sensor_pos(1,j)/1e3 +.2, sensor_pos(2,j)/1e3-.2,lbl_sen,"Color",[0 0 0]);
end
hold on;
% Contour plot of the precision in the plane
[cp,hiso2] = contour(xy_ax,x/1e3,y/1e3,real(GDOPxy)/1e6,contourLevelsxy,'LineColor','k','DisplayName','Precision in the XY-plane');
clabel(cp,hiso2);
legend('Location','NorthEast');
grid off;

% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]'); title('XY GDOP')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    XZ plane    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure();hold on;
plot(sensor_pos(1,:)/1e3,sensor_pos(3,:)/1e3,'o','DisplayName','Sensors','LineWidth',3, ...
    'MarkerFaceColor','blue','MarkerEdgeColor','blue');
for j = 1:nSensors
    lbl_sen = sprintf('S_{%1.0d}',j);
    text(sensor_pos(1,j)/1e3 + .2, sensor_pos(3,j)/1e3 - .2,lbl_sen,"Color",[0 0 0]);
end

% Contour plot of the precision in the plane
[cp,hiso2] = contour(x/1e3,z/1e3,GDOPxz/1e6,contourLevels,'LineColor','k','DisplayName','Precision in the XZ-plane');
clabel(cp,hiso2);
legend('Location','NorthEast');
grid off;
xlim([0 x_span(2)/4e3]); ylim([0 5])
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Altitude [km]'); title('XZ GDOP')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    YZ plane    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure();hold on;
plot(sensor_pos(2,:)/1e3,sensor_pos(3,:)/1e3,'o','DisplayName','Sensors','LineWidth',3, ...
    'MarkerFaceColor','blue','MarkerEdgeColor','blue');
for j = 1:nSensors
    lbl_sen = sprintf('S_{%1.0d}',j);
    text(sensor_pos(2,j)/1e3 +.2, sensor_pos(3,j)/1e3-.2,lbl_sen,"Color",[0 0 0]);
end

% Contour plot of the precision in the plane
[cp,hiso2] = contour(y/1e3,z/1e3,GDOPyz/1e6,contourLevels,'LineColor','k','DisplayName','Precision in the YZ-plane');
clabel(cp,hiso2);
legend('Location','NorthEast');
grid off;
xlim([0 x_span(2)/4e3]);ylim([0 2])

% Adjust the Display
xlabel('Down-range [km]');ylabel('Altitude [km]'); title('YZ GDOP')
hold off


%% Assign output
gdop_xy = GDOPxy;
gdop_xz = GDOPxz;
gdop_yz = GDOPyz;

%% Main sub function for computing plane-wise GDOP

    function [GDOP,axis1,axis2] = gen_GDOPfn(sensor_pos,plane,timingerr,res_a1,res_a2,axis1_min,axis1_max,axis2_min,axis2_max)
    % gdop_fn 
    
    % And computes the Geometric dilution of precision (GDOP) by initializing
    % possible source positions over a volume and computing the ratio of error
    % in position to the error in range for ach pair of sensors.
    
    % DOP represents the positional precision due to the angular separation
    % between the sensors that are used to estimate a source's position.
    
    % Log at end
    % Author: Abeer Nasir Chaudhry®
    
    c = 3e8;
    
    tdoaAccuracy = timingerr;
    
    nSensors = size(sensor_pos,2);
    sensor_pos(3,1:nSensors) = zeros(1,nSensors);
    
    % Create span
    axis1 = linspace(axis1_min,axis1_max,res_a1);
    axis2 = linspace(axis2_min,axis2_max,res_a2);
    
    [X,Y] = meshgrid(axis1,axis2);
    Z = zeros(size(X))+0;
     
    
    % All possible source positions
    targetPos = [X(:) Y(:) Z(:)]';
    
    %% Decentralised sensor pairing
    rxPairs = nchoosek(1:nSensors,2);
    n = size(rxPairs,1);
    
    %% Initialize matrices
    if strcmp(plane,'xy') == 1 || strcmp(plane,'yx') == 1
        invcova11 = zeros(1,numel(X));
        invcova22 = zeros(1,numel(X));
        invcova12 = zeros(1,numel(X));
        GDOP = zeros(size(X));
    else
        invcova11 = zeros(1,numel(Z));
        invcova22 = zeros(1,numel(Z));
        invcova12 = zeros(1,numel(Z));
        GDOP = zeros(size(Z));
    end
    
    
    %% Compute GDOP
    for i = 1:n
        % Access a sensor pair from the network
        s1 = sensor_pos(:,rxPairs(i,1));
        s2 = sensor_pos(:,rxPairs(i,2));
        
        % Calculate directional derivative 
        diff_x1 = targetPos(1,:) - s1(1);
        diff_y1 = targetPos(2,:) - s1(2);
        diff_z1 = targetPos(3,:) - s1(3);
        diff_x2 = targetPos(1,:) - s2(1);
        diff_y2 = targetPos(2,:) - s2(2);
        diff_z2 = targetPos(3,:) - s2(3);
    
        % Range, sensor to target
        r1 = sqrt(diff_x1.^2 + diff_y1.^2 + diff_z1.^2);
        r2 = sqrt(diff_x2.^2 + diff_y2.^2 + diff_z2.^2);
        
        % Directional derivative
        tx = (diff_x1./r1 - diff_x2./r2)./c;
        ty = (diff_y1./r1 - diff_y2./r2)./c;
        tz = (diff_z1./r1 - diff_z2./r2)./c;
    % Compute the inverse of Cov mtx for the 
    invcova11 = invcova11 + tx.^2./tdoaAccuracy;
    invcova22 = invcova22 + ty.^2./tdoaAccuracy;
    invcova12 = invcova12 + tx.*ty./tdoaAccuracy;
    end
    
    % Get the trace of the covariance mtx
    cova11 = invcova22./(invcova11.*invcova22 - invcova12.^2);
    cova22 = invcova11./(invcova11.*invcova22 - invcova12.^2);
    
    
    GDOP(:) = (sqrt(cova11 + cova22));
    
    end
    
    % LOG:
    % 1. Date created (updated): 221110 
end