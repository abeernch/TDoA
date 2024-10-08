%% Ooutput of ELINT SG with 2D Solution and plotting of hyperbola
% The script takes the TOA readings generated by the ELINT SG. The 2D
% solution in xy is determined along with the plotting of resulting
% hyperbolas based on the sensor network and TOA readings.

%% Assumptions and limitations:
% 1. Single emitter
% 2. No multipath effect
% 3. No scan pattern etc implemented. TOA readings based off constant
% amplitude, high SNR IQ data
% 4. Solution extendable to any no of sensors

% The error observed in intersections maybe because of any of the
% following:    
% 1. Height factor
% 2. GDOP
% 3. Error in reported TOA because of low sampling rate

% Next step would be to to the same for xz and yz.

%%
clearvars -except tx rx reported_toa
clc

%% Initialize Sensor Position

sensor_pos = [rx{1}.pos rx{2}.pos rx{3}.pos rx{4}.pos];
nSensors = size(sensor_pos,2);

if nSensors >=3
    % Central Reference Pairing
    tdoaPairs(:,2) = (2:(nSensors))';tdoaPairs(:,1) = 1;

    %% Initialize for solution
    syms x [1 nSensors]
    syms y [1 nSensors]
    syms xx yy c
    syms t [1 nSensors]

    %  Baseline Equations to Solve
    for pairs = 2:size(tdoaPairs,1)
        eqn1        = (sqrt((x(1) - xx)^2 + (y(1) - yy)^2) - sqrt((x(2) - xx)^2 + (y(2) - yy)^2)  == c*(t(1) - t(2)));
        eqn2(pairs) = (sqrt((x(1) - xx)^2 + (y(1) - yy)^2) - sqrt((x(pairs + 1) - xx)^2 + (y(pairs + 1) - yy)^2)  == c*(t(1) - t(pairs + 1)));
    end

    %% Assign parameters to equation variables
    c = 3e8;

    for i=1:nSensors
        temp_varx = strcat( 'x',num2str(i) );
        eval(sprintf('%s = %g',temp_varx,sensor_pos(1,i)));
        temp_vary = strcat( 'y',num2str(i) );
        eval(sprintf('%s = %g',temp_vary,sensor_pos(2,i)));
    end

    % Solve partial equation and create inline functions
    for i = 1:size(tdoaPairs,1) - 1
        sol(i) = solve([eval(eqn1) eval(eqn2(i+1))],[xx yy]);
    end

    sol1Hyp_xx1 = matlabFunction(sol(1).xx(1)); sol1Hyp_xx2 = matlabFunction(sol(1).xx(2));
    sol2Hyp_xx1 = matlabFunction(sol(2).xx(1)); sol2Hyp_xx2 = matlabFunction(sol(2).xx(2));
    sol1Hyp_yy1 = matlabFunction(sol(1).yy(1)); sol1Hyp_yy2 = matlabFunction(sol(1).yy(2));
    sol2Hyp_yy1 = matlabFunction(sol(2).yy(1)); sol2Hyp_yy2 = matlabFunction(sol(2).yy(2));

    compute_gdop = 0; % Enable/disable GDOP calculation
    time_err = 0e-9; % For GDOP analysis

    %% TOA Processing
    for i = 1:nSensors
        obs(i) = length(reported_toa{i});
        min_obs = min(obs);
        if size(reported_toa{i},2) ~= 1
            reported_toa{i} = reported_toa{i}.';
        end
    end
    for i = 1:nSensors
        pri(i,:)     = (1:min_obs-1)./tx{1}.prf;
%         pri(i,:)     = (1:min_obs-1)./diff(reported_toa{i}).';

        toa_arr(i,:) = reported_toa{i}(1:min_obs,1).'./tx{1}.fs;
    end

    toa_arr = [toa_arr(1:nSensors,1) toa_arr(1:nSensors,2:end) - pri];


    %% Initialize Emitter positions
    if tx{1}.platformType == 2
        emitter_pos = tx{1}.pos;
    else
        emitter_pos = tx{1}.trajPos.';
        emitter_pos = emitter_pos(:,1:min_obs);
    end
    nEmitter = size(emitter_pos,2);

    %% Initialize figure and children
    figure;grid on
    %% Plot sensor positions and label
    % Colors
    sens_clr = [0 0.5333 1.0000];
    cm = turbo(nEmitter);
    sens_txtclr = [0 0 0];
    iso_txtclr = [0 0 0];
    em_txtclr = [0 0 0];
    % Adjust Axes
    lim_x = [min([sensor_pos(1,:) emitter_pos(1,:)]) - 10 max([sensor_pos(1,:) emitter_pos(1,:)]) + 10];
    lim_y = [min([sensor_pos(2,:) emitter_pos(2,:)]) - 10  max([sensor_pos(2,:) emitter_pos(2,:)]) + 10];
    xlim(lim_x); ylim(lim_y)
    xlabel('Cross-range (km)');ylabel('Down-range (km)');

    %% Start Processing the TOA readings
    for k = 1:length(toa_arr)
        for n = 1:size(tdoaPairs,1)
            % Draw Isochrone pairs for each reading
            iso{n} = draw_2Disochrone(sensor_pos(:,1),sensor_pos(:,n+1),(toa_arr(1,k) - toa_arr(n+1,k))*c,300e3,20e3);
        end

        %% Find Solution at each reading
        [t1,t2,t3,t4] = deal(toa_arr(1,k),toa_arr(2,k),toa_arr(3,k),toa_arr(4,k));
        %% Evaluate equations
        sol_x(1:4,k) = [sol1Hyp_xx1(t1,t2,t3); sol1Hyp_xx2(t1,t2,t3); sol2Hyp_xx1(t1,t2,t4); sol2Hyp_xx2(t1,t2,t4)];
        sol_y(1:4,k) = [sol1Hyp_yy1(t1,t2,t3); sol1Hyp_yy2(t1,t2,t3); sol2Hyp_yy1(t1,t2,t4); sol2Hyp_yy2(t1,t2,t4)];
        %% Plot isochrones and solution overlays
        for i = 1:nSensors-1
            lbl_iso = sprintf('Isochrone Emitter %1.0d',i);         % Isochrone labels
            iso_obj = plot(iso{i}(1,:),iso{i}(2,:),'LineStyle',':','DisplayName',lbl_iso,'Color',[0.4510 0.9412 0.5490],'LineWidth',1);
            if i ~= 1
                excludeFromLegend(iso_obj);                           % Exclude multiple labels
            end
        end
        % Plot solutions
        for j = 1:4
            sol1 = scatter(sol_x(j,k),sol_y(j,k),'green','diamond','LineWidth',1,'DisplayName','Detections');
        end
        % Plot sensor positions
        for j = 1:(nSensors)
            lbl_sen = sprintf('S_{%1.0d}',j);
            hold on
            sens = scatter(sensor_pos(1,j), sensor_pos(2,j),[],sens_clr,'filled','o','DisplayName','Sensors', ...
                'MarkerEdgeColor','k','LineWidth',0.5);
            text(sensor_pos(1,j) +.2, sensor_pos(2,j)-.2,lbl_sen,"Color",sens_txtclr);
            if j ~=1
                excludeFromLegend(sens);
            end
        end
        ax = gca;
        ax.YDir = 'reverse';
        for j = 1:nSensors - 1
            % TDOA and Sensor lables
            if nEmitter<2
                lbl = sprintf('TDOA_{%1.0d,%1.0d}',j,nSensors);     % Label Isochrones only if single emitter (otherwise clutter)
                txt = text(mean([sensor_pos(1,j),sensor_pos(1,end)]),mean([sensor_pos(2,j),sensor_pos(2,end)])-.2,lbl,"Color",iso_txtclr);
            end
        end
        %% Plot emitter position(s)
        if tx{1}.platformType == 2
            em = scatter(emitter_pos(1),emitter_pos(2),'red','filled','^','DisplayName','Transmitter', ...
                'MarkerEdgeColor','k','LineWidth',0.5);
            lbl_tx = sprintf('T_{%1.0d}',1);
            text(emitter_pos(1,:) +.2,emitter_pos(2,:)+.2,lbl_tx,"Color",em_txtclr);
            if k >1
                excludeFromLegend(em);
            end
        else
            em = scatter(emitter_pos(1,k),emitter_pos(2,k),'red','filled','^','DisplayName','Transmitter', ...
                'MarkerEdgeColor','k','LineWidth',0.5);
            lbl_tx = sprintf('T_{%1.0d}',k);
            text(emitter_pos(1,k) +.2,emitter_pos(2,k)+.2,lbl_tx,"Color",em_txtclr);
            if k >1
                excludeFromLegend(em);
            end
        end


        % Update axes
        drawnow()
%         cla
        sc = findobj('Type','scatter'); txt = findobj('Type','text'); ll = findobj('Type','line');
        delete(sc); delete(txt); delete(ll)
    end
    close
else
    fprintf('At least 03 Sensors required to proceed')
end
% LOG
% 1. Date created 221226 (221229) (230102)