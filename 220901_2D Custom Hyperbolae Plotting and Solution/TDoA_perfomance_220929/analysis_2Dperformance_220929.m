%% Performance analysis of sensor network 
% This script generates TOA readings at a sensor network received from an
% mobile emitter. The TOA readings are used to implement 2D TDOA algorithm
% to get a PF on the mobile emitter. 2D hyperbola are plotted and solutions
% are determined and plotted.

% Author: Abeer Nasir Chaudhry®
% 220929
%%
clear;clc

compute_gdop = 0; % Enable/disable GDOP calculation
time_err = 10e-9; % For GDOP analysis

%% Initialize sensor positions
% sensor_pos = [00e3 -0.20e3 0.20e3  00e3 ;
%               00e3  0.20e3 0.20e3  -0.25e3;
%               00e3 00e3 00e3 00e3];\
% sensor_pos =    [0	3116.05483775860	-11496.8490106353;
%                  0	-7246.75408749039	-21.4179927785758;
%                  0	-75.8919118246145	46.6497395972486];

sensor_pos =    [0	5828.69282906950	3116.05483775860 -11496.8490106353;
                 0	-256.108139194091	-7246.75408749039 -21.4179927785758;
                 0	-12.6655102956660	-75.8919118246145 46.6497395972486];

%% Initialize Source position
[x,y] = circle(sensor_pos(1,1),sensor_pos(2,1),3e3,3);
% emitter_pos = [x;y;ones(1,length(x))*00e3];
emitter_pos  = [100e3;0.5e3;0];
% Emission speed
c = 3e8;

%% TOA Generator function
[nSensors, nEmitter, sensor_pos, emitter_pos, emitter_init, reported_toa, ranges] = mod_ToaGenerator(sensor_pos,emitter_pos,false);

%% Form the Isochrones
% Central sensor pairing (ref sensor = s_1)
tdoaPairs(:,2) = (2:(nSensors))';
tdoaPairs(:,1) = 1;

for n = 1:nEmitter
    for i = 1:size(tdoaPairs,1)
            iso{n,i} = draw_2Disochrone(sensor_pos(:,1),sensor_pos(:,i+1),(ranges(1,n)-ranges(i+1,n)),100e3,25e3);
    end
end

%% Compute GDOP
if compute_gdop == 1
    [gdop,xspan,yspan] = gdop_fn(sensor_pos,time_err);
end
%% Find the solution to the hyperbolic eqns for 2D TDoA
for n = 1:nEmitter
    for i = 1:nSensors % Assign positions to new variables
        x_s(i) = sensor_pos(1,i);
        y_s(i) = sensor_pos(2,i);
        t(i) = reported_toa(i,n)% + (rand(1,1)*(100e-9)) - (15e-9);
    end

    %% Evaluate equations
    for pairs = 2:size(tdoaPairs,1)
        syms x y
        s1 = eval(sqrt((x_s(1) - x)^2 + (y_s(1) - y)^2) - sqrt((x_s(2) - x)^2 + (y_s(2) - y)^2)  == c*(t(1) - t(2)));
        s2 = eval(sqrt((x_s(1) - x)^2 + (y_s(1) - y)^2) - sqrt((x_s(pairs + 1) - x)^2 + (y_s(pairs + 1) - y)^2)  == c*(t(1) - t(pairs + 1)));
%         s2 = eval(sqrt((x_s(1) - x)^2 + (y_s(1) - y)^2) - sqrt((x_s(3) - x)^2 + (y_s(3) - y)^2)  == c*(t(1) - t(3)));
        
        % Solve eqns for (x,y) given the values of other eqn variables 
        sol = solve([s1 s2],[x y]);
%         xsol(1:2,:) = sol.x;
%         ysol(1:2,:) = sol.y;
        xsol(1:2,pairs - 1) = sol.x;
        ysol(1:2,pairs - 1) = sol.y;
    end
    
    % Compile solution
    solu{n} = [(eval(xsol(:))).';(eval(ysol(:))).'];
    
    %% Reject ghost detections.
    % Extract correct solution on the basis of frequency of occurence from
    % different baselines. 
    final_sol{n} = round(solu{n},2);
%     final_sol{n} = final_sol{n}(final_sol{n}~=0);
%     final_sol{n} = reshape(final_sol{n},2,length(final_sol{n})/2);
    final_sol{n} = mode(final_sol{n},2);

    % Print results
    detect = sprintf('Detection %1.0d:\n x (km): %2.5f\n y (km): %2.5f\n', n,final_sol{n}(1,:)./1e3,final_sol{n}(2,:)./1e3);
end

%% Plotting the results and the solution(s)
figure;grid on
if compute_gdop == 1 
%% Plotting GDOP
    img = imagesc(gdop./1e6,'XData',xspan./1e3,'YData',yspan./1e3);
    colorbar
    set(gca,'YDir','normal')
    sens_clr = [1 1 1];
    % Apply colour jet for multiple emitters
    cm = jet(nEmitter);
    sens_txtclr = [1 1 1];
    iso_txtclr = [0.9412 0.9412 0.9412];
    em_txtclr = [1 0 0];
else
    sens_clr = [0 0.5333 1.0000];
    cm = turbo(nEmitter);
    sens_txtclr = [0 0 0];
    iso_txtclr = [0 0 0];
    em_txtclr = [0 0 0];
end
%% Plot sensor positions and label
% Plot sensor positions
hold on
for j = 1:(nSensors)
    lbl_sen = sprintf('S_{%1.0d}',j);
    hold on
    % Plot sensor positions
    sens = scatter(sensor_pos(1,j)/1e3, sensor_pos(2,j)/1e3,[],sens_clr,'filled','o','DisplayName','Sensors', ...
                   'MarkerEdgeColor','k','LineWidth',0.5);
    text(sensor_pos(1,j)/1e3 +.2, sensor_pos(2,j)/1e3-.2,lbl_sen,"Color",sens_txtclr);
    if j ~=1
        excludeFromLegend(sens);
    end
end
%% Plot emitter position(s)
hold on
for i = 1:nEmitter
    em = scatter(emitter_pos(1,i)/1e3,emitter_pos(2,i)/1e3,'red','filled','^','DisplayName','Transmitter', ...
                 'MarkerEdgeColor','k','LineWidth',0.5);
    lbl_tx = sprintf('T_{%1.0f}',i);
    text(emitter_pos(1,i)/1e3 +.2,emitter_pos(2,i)/1e3+.2,lbl_tx,"Color",em_txtclr);
    if i ~=1
        excludeFromLegend(em);
    end
end

%% Plot isochrones and solution overlays
for n = 1:nEmitter
    for i = 1:nSensors-1
        lbl_iso = sprintf('Isochrone Emitter %1.0d',n);         % Isochrone labels
        iso_obj = plot(iso{n,i}(1,:)/1e3,iso{n,i}(2,:)/1e3,'LineStyle',':','DisplayName',lbl_iso,Color=cm(n,:));
        if i ~=1
            excludeFromLegend(iso_obj);                           % Exclude multiple labels
        end
        % TDOA and Sensor lables
        if nEmitter<2
            lbl = sprintf('TDOA_{%1.0d,%1.0d}',i,nSensors);     % Label Isochrones only if single emitter (otherwise clutter)
            text(mean([sensor_pos(1,i)/1e3,sensor_pos(1,end)/1e3]),mean([sensor_pos(2,i)/1e3,sensor_pos(2,end)/1e3])-.2,lbl,"Color",iso_txtclr);
        end
    end
    % Plot solutions
    sol_obj = scatter(solu{n}(1,:)/1e3,solu{n}(2,:)/1e3,'green','diamond','LineWidth',1,'DisplayName','Detections');
    if n > 1
        excludeFromLegend(sol_obj);                           % Exclude multiple labels
    end

end

hold on
% Adjust Axes
lim_x = [min([sensor_pos(1,:)/1e3 emitter_pos(1,:)/1e3]) - 10 max([sensor_pos(1,:)/1e3 emitter_pos(1,:)/1e3]) + 10];
lim_y = [min([sensor_pos(2,:)/1e3 emitter_pos(2,:)/1e3]) - 10  max([sensor_pos(2,:)/1e3 emitter_pos(2,:)/1e3]) + 10];
xlim(lim_x); ylim(lim_y)

xlabel('Cross-range (km)');ylabel('Down-range (km)');
% legend('Location','northeast');


%% False solutions identifier
for i = 1:nEmitter
    sol_chk = real(round(solu{i},2));
    chk = sol_chk(:,1) ~= sol_chk(:,:);
    if any(any(chk))
        error_pos(i) = i;
    end
end

% LOG
% 1. Date created (updated): 220927 (221003) (221004)
% 2. GDOP function added
% 3. Dotted lines changed to solid lines