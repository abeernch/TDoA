function [invGDOP] = gdop3_fn(sensor_pos,timingerr,grd_res,x_min,x_max,y_min,y_max)
% gdop_fn 
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
    sensor_pos = [00e3 -25e3 25e3 00e3;00e3 25e3 25e3 -25e3; 1e3 -0.50e3 0.150e3 00e3];
%     sensor_pos = [00e3 -25e3 25e3 00e3;00e3 25e3 25e3 -25e3; 00e3 00e3 00e3 00e3];
% sensor_pos = [00e3 00e3 00e3 00e3;00e3 25e3 25e3 -25e3; 1e3 -0.50e3 0.150e3 00e3];

    timingerr = 30e-9;
    grd_res = 1e2;
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 2
    timingerr = 30e-9;
    grd_res = 1e2;
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;
    

elseif nargin < 3
    grd_res = 1e2;
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;
    

elseif nargin < 4
    x_min = -200e3;
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;
    

elseif nargin < 5
    x_max = 200e3;
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 6
    y_min = -200e3;
    y_max = 200e3;

elseif nargin < 7
    y_max = 200e3;
end


c = 3e8;
tdoaAccuracy = timingerr;
nSensors = size(sensor_pos,2);

nDim = size(sensor_pos,1);
if nDim ~= 3
    error(['Enter sensor positions in the following format [x1 x2 x3...xn; ' ...
                                                          ' y1 y2 y3...yn; ' ...
                                                          ' z1 z2 z3...zn])'])
end

% Create expanse
x = linspace(x_min,x_max,grd_res);
y = linspace(y_min,y_max,grd_res);
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

% All possible source positions
targetPos = [X(:) Y(:) Z(:)]';

%% Decentralised sensor pairing
rxPairs = nchoosek(1:nSensors,2);
n = size(rxPairs,1);

%% Initialize matrices
invcovxx = zeros(1,numel(X));
invcovyy = zeros(1,numel(X));
invcovxy = zeros(1,numel(X));
invGDOP = zeros(size(X));

%% Compute GDOP
for i = 1:n
    dr1 = sensor_pos(:,rxPairs(i,1));
    dr2 = sensor_pos(:,rxPairs(i,2));
    dx1 = targetPos(1,:) - dr1(1);
    dy1 = targetPos(2,:) - dr1(2);
    dz1 = targetPos(3,:) - dr1(3);
    dx2 = targetPos(1,:) - dr2(1);
    dy2 = targetPos(2,:) - dr2(2);
    dz2 = targetPos(3,:) - dr2(3);

    r1 = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
    r2 = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
    
    tx = (dx1./r1 - dx2./r2)./c;
    ty = (dy1./r1 - dy2./r2)./c;
    
    invcovxx = invcovxx + tx.^2./tdoaAccuracy;
    invcovyy = invcovyy + ty.^2./tdoaAccuracy;
    invcovxy = invcovxy + tx.*ty./tdoaAccuracy;
end

covxx = invcovyy./(invcovxx.*invcovyy - invcovxy.^2);
covyy = invcovxx./(invcovxx.*invcovyy - invcovxy.^2);
invGDOP(:) = (sqrt(covxx + covyy));

%% Plotting
figure;
%% Plotting GDOP
    img = imagesc(gdop_xy./1e6,'XData',x./1e3,'YData',y./1e3);
    grid on; xlabel('Cross-range (km)'); ylabel('Down-range (km)'); title('GDOP');
    colorbar
    set(gca,'YDir','normal')
    sens_clr = [1 1 1];
    % Apply colour jet for multiple emitters
    sens_txtclr = [1 1 1];
    iso_txtclr = [0.9412 0.9412 0.9412];
    em_txtclr = [1 0 0];

%% Plot sensor positions
hold on
for j = 1:(nSensors)
    lbl_sen = sprintf('S_{%1.0d}',j);
    hold on
    % Plot sensor positions
    sens = scatter3(sensor_pos(1,j)/1e3, sensor_pos(2,j)/1e3, sensor_pos(3,j)/1e3,[],sens_clr,'filled','o','DisplayName','Sensors', ...
                   'MarkerEdgeColor','k','LineWidth',0.5);
    text(sensor_pos(1,j)/1e3 +.2, sensor_pos(2,j)/1e3-.2, sensor_pos(3,j)/1e3,lbl_sen,"Color",sens_txtclr);
    if j ~=1
        excludeFromLegend(sens);
    end
end


end


% LOG:
% 1. Date created (updated): 221004 (221105)
