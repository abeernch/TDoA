%% Least Squares Augmented 2D Hyperbolic Localization using TDoA measurements from ADS-B

% Author: Abeer Nasir Chaudhry®
% 250131

%% Initialize
clear;close all; clc
%% Configure the simulation
% Select file to load
filename = ('log_4074e4_complete_4nodes_31mar_partial_working.mat');
load(filename);

% Annotate errors on plot? 
annotateErr = 0; steps = 10;

% Whether to plot isochrones and platforms positions or just get the
% solutions
plotDo = 1;
% Whether to plot on local axes or geoaxes (with geo coords)
localPlot = 1;
analyzeData = 0;


%% Convert Sensor LAT,LON to local
% Node Positions in geographic coordinates
% Node 1: NASTP (reference sensor)
% Node 2: Swana
% Node 3: I-14
% Node 4: ChakShahzad
refSen = 1;
[lat_node{1},lon_node{1},alt_node{1}] = deal(33.609144,73.104425,510);
[lat_node{2},lon_node{2},alt_node{2}] = deal(33.5438077,73.1379735,439);
[lat_node{3},lon_node{3},alt_node{3}] = deal(33.6088889,72.9805555,567);
[lat_node{4},lon_node{4},alt_node{4}] = deal(33.674690,73.184789,536);

nSensors = size(lat_node,2);

% Initialize Sensor Position
% Pre-allocate
rx = cell(1,nSensors);

for i = 1:nSensors
    [rx{i}.pos(1,1),rx{i}.pos(2,1),rx{i}.pos(3,1)] = latlon2local(lat_node{i},lon_node{i},alt_node{i},[lat_node{1},lon_node{1},alt_node{1}]);
end

sensor_pos = [rx{1}.pos (rx{2}.pos) rx{3}.pos rx{4}.pos];

c = physconst('Lightspeed');

%% Call Cleaner Function for raw recorded observations
[toa_arr(1,:),toa_arr(2,:),toa_arr(3,:),toa_arr(4,:)] = deal(TOA1,TOA2,TOA3,TOA4);

% Outliers in the TDoA are removed by computing the MAD of the entire data
% about a window median. If the sample lies outside the 3sigma limit (3
% times the std dev of the data in the window), it is replaced by the
% median of the window.
[latSourceClean,lonSourceClean,altSourceClean,tdoaClean] = fn_tdoacleaner3(nSensors,toa_arr,lat1,long1,alt1,analyzeData);

%% Initialize Emitter positions
tx{1}.platformType = 1;
[txpos_x,txpos_y,txpos_z] = latlon2local(latSourceClean,lonSourceClean,altSourceClean/3.281,[lat_node{1},lon_node{1},alt_node{1}]);
tx{1}.trajPos = [txpos_x.';txpos_y.';txpos_z.'];
emitter_pos = tx{1}.trajPos;
nEmitter = size(emitter_pos,2);

%% Central Reference System Pairing
tdoaPairs(:,2) = (2:(nSensors))';tdoaPairs(:,1) = 1;

%% Initialize for solution
syms x1 x2 x3 x4 xx y1 y2 y3 y4 yy z1 z2 z3 z4 zz c tdoa12 tdoa13 tdoa14
eqn1 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x2-xx)^2 + (y2-yy)^2 + (z2-zz)^2)  == c*(tdoa12);
eqn2 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x3-xx)^2 + (y3-yy)^2 + (z3-zz)^2)  == c*(tdoa13);
eqn3 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x4-xx)^2 + (y4-yy)^2 + (z4-zz)^2)  == c*(tdoa14);

[x1,y1,z1]       = deal(sensor_pos(1,1),sensor_pos(2,1),sensor_pos(3,1));
[x2,y2,z2]       = deal(sensor_pos(1,2),sensor_pos(2,2),sensor_pos(3,2));
[x3,y3,z3]       = deal(sensor_pos(1,3),sensor_pos(2,3),sensor_pos(3,3));
[x4,y4,z4]       = deal(sensor_pos(1,4),sensor_pos(2,4),sensor_pos(3,4));

c = physconst('Lightspeed');

% Solve partial equation and create inline functions
sol = solve([eval(eqn1) eval(eqn2) eval(eqn3)],[xx yy zz]);

solHyp_x1 = matlabFunction(sol.xx(1));solHyp_x2 = matlabFunction(sol.xx(2));
solHyp_y1 = matlabFunction(sol.yy(1));solHyp_y2 = matlabFunction(sol.yy(2));
solHyp_z1 = matlabFunction(sol.zz(1));solHyp_z2 = matlabFunction(sol.zz(2));

%=======================================================================================================%
%% Begin Plotting and updating
if plotDo == 1
    if localPlot == 1
        figure;grid on
        ax = gca;
        ax.Title.String = '04 Node TDoA Based Hyperbola Position Fix';
        ax.XLabel.String = 'Cross-range (m)';
        ax.YLabel.String = 'Down-range (m)';
        ax.XDir = 'reverse';
        xlim([-200e3 200e3]); ylim([-200e3 200e3])
        legend

        %% Plot sensor positions and label Colors
        sens_clr = [0 0.5333 1.0000];
        cm = turbo(nEmitter);
        sens_txtclr = [0 0 0];
        iso_txtclr = [0 0 0];
        em_txtclr = [0 0 0];

        for j = 1:(nSensors)
            lbl_sen = sprintf('S_{%1.0d}',j);
            hold on
            sens = scatter(sensor_pos(1,j), sensor_pos(2,j),100,sens_clr,'filled','o','DisplayName','Sensors', ...
                'MarkerEdgeColor','k','LineWidth',1);
            text(sensor_pos(1,j) + 500, sensor_pos(2,j) + 1000,lbl_sen,"Color",sens_txtclr,'FontWeight','bold');
            if j ~=1
                excludeFromLegend(sens);
            end
        end

        % Pre-allocate solution arrays
        sol_x  = zeros(2,length(tdoaClean{1}));
        sol_y  = zeros(2,length(tdoaClean{1}));

        % Initialize for logging and appending data from previous batch
        acc_sol1X = [];acc_sol1Y = [];
        acc_sol2X = [];acc_sol2Y = [];
        acc_lse = [];acc_lse = [];
        acc_trailX = [];acc_trailY = [];
        inc = 1;

        for k = 1:length(tdoaClean{1})
            %% Evaluate equations
            sol_x(1:2,k) = [solHyp_x1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_x2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
            sol_y(1:2,k) = [solHyp_y1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_y2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
            sol_z(1:2,k) = [solHyp_z1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_z2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];            
            pos_est = [sol_x(2,k).' sol_y(2,k).'];
            
            measuredRDoA = [tdoaClean{1}(k) tdoaClean{2}(k) tdoaClean{3}(k)].'.*c;
            %% Implement the OLS solution
            [error, posPredict] = fn_tdoaOLS(sensor_pos, refSen, pos_est.',measuredRDoA, 500e-9, 1e-12, 50);
    
            % LSest(:,k) = posPredict(:,end);
            % err(k) = error(end);

            %% Compute 10-point RMS error in the real-time
            if annotateErr == 1
                if mod(k,steps) == 0
                    realTerr = vecnorm(emitter_pos(1:2,inc : inc + steps - 1) - [sol_x(1,inc : inc + steps - 1);sol_y(1,inc : inc + steps - 1)]);
                    rmsE = sqrt(sum(abs(realTerr.^2))/length(realTerr));
                    err_mkr = text(emitter_pos(1,inc(end))-5e2,emitter_pos(2,inc(end))-1.5e3,string(rmsE) + ' m','FontSize',6,'Color',[0.45 0.45 0.45]);
                    inc = inc + steps;
                end
            end
                      
            %% Plot emitter position(s)
            if k <= 1
                em = scatter(emitter_pos(1,k),emitter_pos(2,k),100,'red','filled','^','DisplayName','Transmitter', ...
                    'MarkerEdgeColor','k','LineWidth',1.5);
                trail = scatter(emitter_pos(1,k),emitter_pos(2,k),2,'black','filled','^','DisplayName','Emitter Location History','MarkerEdgeColor','k','LineWidth',0.1);
            else
                em.XData = emitter_pos(1,k);
                em.YData = emitter_pos(2,k);
                trail.XData = [acc_trailX emitter_pos(1,k)];
                trail.YData = [acc_trailY emitter_pos(2,k)];
                acc_trailX = trail.XData;
                acc_trailY = trail.YData;
            end

            if k >1
                excludeFromLegend(trail);
            end
            lbl_tx = sprintf('T_{%1.0d}',k);
            tx_t = text(emitter_pos(1,k) + 500,emitter_pos(2,k)+1000,lbl_tx,"Color",em_txtclr);

            %% Plot solution overlays
            % Solution overlays
            if k <= 1
                sol1 = plot(sol_x(1,k),sol_y(1,k),'LineWidth',0.5,'Color',[0.7 0.1 0.1 ],'LineStyle','none','Marker','.','MarkerSize',15,'DisplayName','TDoA Position Fix');
                % sol2 = plot(sol_x(2,k),sol_y(2,k),'LineWidth',0.5,'Color',[1.0000 0.4118 0.1608],'LineStyle','none','Marker','.','MarkerSize',15,'DisplayName','TDoA Position Fix');
                est = scatter(posPredict(1,end), posPredict(2,end),'MarkerFaceColor','g','Marker','+','DisplayName','LSE');
            else
               est.XData = [acc_lse real(posPredict(1,end))]; est.YData = [acc_lse real(posPredict(2,end))]; 
                % sol2.XData = [acc_sol2X real(sol_x(1,k))];
                % sol2.YData = [acc_sol2Y real(sol_y(1,k))];
                sol1.XData = [acc_sol1X real(sol_x(2,k))];
                sol1.YData = [acc_sol1Y real(sol_y(2,k))];
                % acc_sol2X = real(sol2.XData); acc_sol2Y = real(sol2.YData);
                acc_sol1X = real(sol1.XData); acc_sol1Y = real(sol1.YData);
            end
          
            %% Update axes
            drawnow()
            delete(tx_t);
        end
    else
        map = importdata('map_expanded.png');
        %% Map

        % Maphighres
        minlon = 70.37154;
        maxlon = 77.25746;
        minlat = 29.53471;
        maxlat = 35.52732;
        scalex = 8192/(maxlon - minlon);
        scaley = 8192/(maxlat - minlat);

        % Pre-allocate Geo-Sensor position variable
        Geosensor_pos = zeros(3,nSensors);

        for i = 1:nSensors
            posx = (lon_node{i} - minlon)*scalex;
            posy = abs((lat_node{i} - maxlat)*scaley);
            Geosensor_pos(1:3,i) =[posx; posy; alt_node{i}];
        end

        posx = (long1 - minlon)*scalex;
        posy = abs((lat1 - maxlat)*scaley);
        Geoemitter_pos = [posx.';posy.'];

        figure;
        gx = gca;
        imshow(map.cdata);
        gx.Title.String = '04 Node TDoA Based Hyperbola Position Fix';
        legend

        %% Plot sensor positions and label Colors
        sens_clr = [0 0.5333 1.0000];
        cm = turbo(nEmitter);
        sens_txtclr = [0 0 0];
        iso_txtclr = [0 0 0];
        em_txtclr = [0 0 0];

        for j = 1:(nSensors)
            lbl_sen = sprintf('S_{%1.0d}',j);
            hold on
            sens = scatter(Geosensor_pos(1,j), Geosensor_pos(2,j),100,sens_clr,'filled','o','DisplayName','Sensors', ...
                'MarkerEdgeColor','k','LineWidth',1.5);
            text(Geosensor_pos(1,j) + 15, Geosensor_pos(2,j)+ 15,lbl_sen,"Color",sens_txtclr,"FontWeight","bold");
            if j ~=1
                excludeFromLegend(sens);
            end
        end
        hold(gx,"on")
        axis equal

        % Pre-allocate isochrone, solution
        iso = cell(1,size(tdoaPairs,1));
        Geoiso = cell(1,size(tdoaPairs,1));
        sol_x  = zeros(2,length(tdoaClean{1}));
        sol_y  = zeros(2,length(tdoaClean{1}));
        lbl_iso = sprintf('Isochrone');
        acc_sol1X = [];acc_sol1Y = [];
        acc_sol2X = [];acc_sol2Y = [];
        acc_trailX = [];acc_trailY = [];
        inc = 1;
        %% begin processing loop
        for k = 1:length(tdoaClean{1})
            for n = 1:size(tdoaPairs,1)
                % Draw Isochrone pairs for each reading
                iso{n} = draw_2Disochrone(sensor_pos(:,1),sensor_pos(:,n+1),tdoaClean{n}(k)*c,100e3,50e3);
                [isoy,isox] = local2latlon(real(iso{n}(1,:)),real(iso{n}(2,:)),0,[lat_node{1},lon_node{1},0]);
                isoxx = (isox - minlon)*scalex;
                isoyy = abs((isoy - maxlat)*scaley);
                Geoiso{n} = [isoxx;isoyy];
            end

            %% Find Solution at each reading
            %% Evaluate equations
            sol_x(1:2,k) = [solHyp_x1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_x2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
            sol_y(1:2,k) = [solHyp_y1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_y2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
            sol_z(1:2,k) = [solHyp_z1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_z2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
            [solx_1,soly_1,~] = local2latlon(real(sol_x(1,k)),real(sol_y(1,k)),0,[lat_node{1},lon_node{1},0]);
            [solx_2,soly_2,~] = local2latlon(real(sol_x(2,k)),real(sol_y(2,k)),0,[lat_node{1},lon_node{1},0]);
            solxx1 = (soly_1 - minlon)*scalex;
            solxx2 = (soly_2 - minlon)*scalex;
            solyy1 = abs((solx_1 - maxlat)*scaley);
            solyy2 = abs((solx_2 - maxlat)*scaley);
            Geosol_x(1:2,k) = [solxx1;solxx2];
            Geosol_y(1:2,k) = [solyy1;solyy2];
            pos_est = [Geosol_x(2,k).' Geosol_y(2,k).'];

            %% Compute 10-point RMS error in the real-time
            if annotateErr == 1
                if mod(k,steps) == 0
                    realTerr = vecnorm(emitter_pos(1:2,inc : inc + steps - 1) - [sol_x(1,inc : inc + steps - 1);sol_y(1,inc : inc + steps - 1)]);
                    rmsE = sqrt(sum(abs(realTerr.^2))/length(realTerr));
                    err_mkr = text(Geoemitter_pos(1,inc(end))+10,Geoemitter_pos(2,inc(end))+5,string(rmsE) + ' m','FontSize',6,'Color',[0.45 0.45 0.45]);
                    inc = inc + steps;
                end
            end
            %% Plot emitter position(s)
            if k <= 1
                em = scatter(Geoemitter_pos(1,k),Geoemitter_pos(2,k),100,'red','filled','^','DisplayName','Transmitter', ...
                    'MarkerEdgeColor','k','LineWidth',1.5);
                trail = scatter(Geoemitter_pos(1,k),Geoemitter_pos(2,k),2,'black','filled','^','DisplayName','Emitter Location History','MarkerEdgeColor','k','LineWidth',0.1);
            else
                em.XData = Geoemitter_pos(1,k);
                em.YData = Geoemitter_pos(2,k);
                trail.XData = [acc_trailX Geoemitter_pos(1,k)];
                trail.YData = [acc_trailY Geoemitter_pos(2,k)];
                acc_trailX = trail.XData;
                acc_trailY = trail.YData;
            end
            lbl_tx = sprintf('T_{%1.0d}',k);
            tx_t = text(Geoemitter_pos(1,k)+15,Geoemitter_pos(2,k)+15,lbl_tx,"Color",em_txtclr,"FontWeight","bold");

            %% Plot isochrones and solution overlays
            if k <= 1
                sol1 = plot(Geosol_x(1,k),Geosol_y(1,k),'LineWidth',1,'Color',[0.7 0.1 0.1],'LineStyle','none','Marker','.','MarkerSize',15,'DisplayName','TDoA Position Fix');
                sol2 = plot(Geosol_x(2,k),Geosol_y(2,k),'LineWidth',1,'Color',[1.0000 0.4118 0.1608],'LineStyle','none','Marker','.','MarkerSize',15,'DisplayName','TDoA Position Fix');
            else
                sol2.XData = [acc_sol2X real(Geosol_x(1,k))];
                sol2.YData = [acc_sol2Y real(Geosol_y(1,k))];
                sol1.XData = [acc_sol1X real(Geosol_x(2,k))];
                sol1.YData = [acc_sol1Y real(Geosol_y(2,k))];
                acc_sol2X = real(sol2.XData); acc_sol2Y = real(sol2.YData);
                acc_sol1X = real(sol1.XData); acc_sol1Y = real(sol1.YData);
            end

            for i = 1:nSensors-1
                if k <= 1
                    iso_obj(i) = plot((Geoiso{i}(1,:)),(Geoiso{i}(2,:)),'LineStyle','-.','DisplayName',lbl_iso,'Color',[0.65 0.65 0.65],'LineWidth',1.5);
                else
                    iso_obj(i).XData = Geoiso{i}(1,:);
                    iso_obj(i).YData = Geoiso{i}(2,:);
                end

                if k > 1
                    excludeFromLegend(iso_obj(1));                           % Exclude multiple labels
                    excludeFromLegend(iso_obj(2));
                end

            end
            if trackDo == 1
                if k == 1
                    track_plot = plot(real(trk_pos(1)),real(trk_pos(2)),'LineWidth',1,'Color',[0.3922 0.8314 0.0745],'DisplayName','Track History','Marker','.','LineStyle','none');
                end
                track_plot.XData = ([acc_trkX trk_pos(1,:)]);
                track_plot.YData = ([acc_trkY trk_pos(2,:)]);
                acc_trkX = track_plot.XData;
                acc_trkY = track_plot.YData;
            end
            %% Update axes
            drawnow()
            delete(tx_t);
        end
    end
else
    %% Only determine the solution (No plotting and visuals)
    % Pre-allocate isochrone, solution
    sol_x  = zeros(2,length(tdoaClean{1}));
    sol_y  = zeros(2,length(tdoaClean{1}));

    for k = 1:length(tdoaClean{1})
        %% Find Solution at each reading
        %% Evaluate equations
        sol_x(1:2,k) = [solHyp_x1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_x2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
        sol_y(1:2,k) = [solHyp_y1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_y2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
        sol_z(1:2,k) = [solHyp_z1(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k)); solHyp_z2(tdoaClean{1}(k),tdoaClean{2}(k),tdoaClean{3}(k))];
        pos_est = [real(sol_x(2,k).') real(sol_y(2,k).')];

     end
    % plot computed stuff
    figure;scatter(sensor_pos(1,:),sensor_pos(2,:),[100],'filled')
    hold on
    scatter(emitter_pos(1,:),emitter_pos(2,:),[15],'filled','k')
    plot(sol_x(2,:),sol_y(2,:),'.','Color',[1.0000 0.4118 0.1608])
    plot(sol_x(1,:),sol_y(1,:),'.','Color',[0.7 0.1 0.1])
end

%% Compute Error in 2D Position Fix
%% Solution 1
sol1 = [sol_x(1,:); sol_y(1,:)];
err_x1 = abs(emitter_pos(1,:) - sol1(1,:));
err_y1 = abs(emitter_pos(2,:) - sol1(2,:));
figure;
subplot(3,1,1); scatter(1:length(err_x1),vecnorm(emitter_pos(1:2,:) - sol1),15,'filled'); title('Absolute Range Error'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(3,1,2); scatter(1:length(err_x1),err_x1,15,'filled'); title('Error in X'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(3,1,3); scatter(1:length(err_y1),err_y1,15,'filled'); title('Error in Y'); xlabel('TOA Sample #'); ylabel('Error (meters)')
sgtitle('Error: Solution 1')

%% Solution 2
sol2 = [sol_x(2,:); sol_y(2,:)];
err_x2 = abs(emitter_pos(1,:) - sol2(1,:));
err_y2 = abs(emitter_pos(2,:) - sol2(2,:));
figure;
subplot(3,1,1); scatter(1:length(err_x2),vecnorm(emitter_pos(1:2,:) - sol2),15,'filled'); title('Absolute Range Error'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(3,1,2); scatter(1:length(err_x2),err_x2,15,'filled'); title('Error in X'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(3,1,3); scatter(1:length(err_y2),err_y2,15,'filled'); title('Error in Y'); xlabel('TOA Sample #'); ylabel('Error (meters)')
sgtitle('Error: Solution 2')

%% Compute TDoA RMSE for both solutions (Ommited outliers)
absErr1 = vecnorm(emitter_pos(1:2,:) - sol1);
absErr2 = vecnorm(emitter_pos(1:2,:) - sol2);
absErr1_new = absErr1(absErr1<100e3);
absErr2_new = absErr2(absErr2<100e3);


rmse1 = sqrt(sum(abs(absErr1_new.^2))/length(absErr1_new))
rmse2 = sqrt(sum(abs(absErr2_new.^2))/length(absErr2_new))
r_max = vecnorm(tx{1}.trajPos(1:2,1) - sensor_pos(1:2,1));
r_min = vecnorm(tx{1}.trajPos(1:2,end) - sensor_pos(1:2,1));
perc_err1 = (rmse1/((r_max+r_min)/2))*100
perc_err2 = (rmse2/((r_max+r_min)/2))*100

%% Compute Tracker RMSE
if localPlot == 0
    sol1 = [Geosol_x(1,:); Geosol_y(1,:)];
    sol2 = [Geosol_x(2,:); Geosol_y(2,:)];
    if trackDo == 1
        TrackErr = vecnorm(sol2 - trk_pos);
        TrackRmse = sqrt(sum(abs(TrackErr.^2))/length(TrackErr))

        figure;subplot(3,1,1);plot(1:length(TrackErr),TrackErr,'.','MarkerSize',10)
        title('Tracker Absolute Error')
        ylabel('Error (m)')
        xlabel('TOA Sample #')
    end
else
end
%% LOG:
% Created: 230206 