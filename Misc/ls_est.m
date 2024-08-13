% Least-Squares Solution for Determination of Emitter Location using TDOA
clear;clc;close all

%% Universal Vars
c = physconst('Lightspeed');

%% Initialize Sensor Network Node Positions

grid on;
axis([-200 200 -200 200])
[xs,ys] = ginput();
nSensors = length(xs)+1;
close;

% pos_sensor = [xs.' ; ys.'];
pos_sensor = cat(2,[xs.' ; ys.'],zeros(2,1))*1e3;
%% Define Sensor Performance
timingError = 1e-9;                 % Sensor Timing Error
rangeError = timingError*c;         % Corresponding range error

Ctoa = timingError^2*eye(nSensors); % TOA Error Covariance matrix
Croa = rangeError^2*eye(nSensors);  % Range Error Covariance matrix

%% Initialize Emitter Location
grid on;
axis([-200 200 -200 200])
[xe,ye] = ginput(1);
pos_emitter = [xe.'; ye.']*1e3;
close;

%% Calculate Difference in Ranges of Arrival among Sensors
roa = range_measurement(pos_sensor,pos_emitter);

%% Generate Range Measurement Errors
montecarlo_iters = 1e2;

rng_offset = rangeError*randn(nSensors,montecarlo_iters);      % [m]
diffRange_error = rng_offset(1:end-1,:) - rng_offset(end,:);   % Generate differential range msmnt noise
noisy_range = roa + diffRange_error;                           % Noisy range 

%% Initial estimate of the emitters position
grid on;
axis([-200 200 -200 200])
[x_init,y_init] = ginput(1); % km
emitter_init = [x_init.'; y_init.']*1e3;
close;
%% Least-Squares Estimation Algorithm
% Initialize variables
iters = 500;
epsilon = 100;                  % desired search resolution size
est = zeros(2,iters);           %
current_est = zeros(2,iters);   %
cov_mtx = zeros(2,2,iters);     % Error Cov Mtx for 

wb = waitbar(0,'Performing Monte Carlo Simulation...');
for i = 1 : montecarlo_iters
    waitbar(i/montecarlo_iters,wb,'Performing Monte Carlo Simulation...')

    % Compute solutions
    [~,est_iter] = compute_least_square(pos_sensor,noisy_range(:,i),Croa,emitter_init,epsilon,iters,true,false,[]);
    
    num_iter_ls = size(est_iter,2);

    if num_iter_ls > iters
        est_iter = est_iter(:,1:iters);
        num_iter_ls = iters;
    end

    current_est(:,1:num_iter_ls) = est_iter;

    if num_iter_ls < iters
        current_est(:,num_iter_ls+1:iters) = est_iter(:,end);
    end

    if i==1
        x_ls = current_est;
    end

    % Update Error
    error_update = pos_emitter - current_est;

    % Compute updated error Cov Mtx
    cov_mtx = cov_mtx + reshape(error_update,2,1,iters).*reshape(error_update,1,2,iters)/montecarlo_iters;
end
fprintf('Complete.\n');

F = findall(0,'type','figure');
delete(F)

plotIndex = [1:10,20:20:100,200:200:iters];
est_plot = x_ls(:,plotIndex,1);

%% Compute CEP50
cep50 = zeros(1,iters);
for ii = 1 : iters
    cep50(ii) = computeCEP50(cov_mtx(:,:,ii))/1e3; % [km]
end
%% Compute CRLB
crlb = computeCRLB(pos_sensor,pos_emitter,Ctoa);
crlb_cep50 = computeCEP50(crlb)/1e3; % [km]
crlb_ellipse = drawErrorEllipse(pos_emitter,crlb,100,90);

%% Plotting Results
sub_sampling = min(20,numel(plotIndex)); % Sub-sampling to plot results

fig_geo_a = figure;
plot(pos_sensor(1,:)/1e3,pos_sensor(2,:)/1e3,'ko','DisplayName','Sensors','LineWidth',1);
hold on;
plot(pos_emitter(1,:)/1e3,pos_emitter(2,:)/1e3,'k^','DisplayName','Emitter','LineWidth',1);
plot(est_plot(1,1:floor(sub_sampling/2),1)/1e3,est_plot(2,floor(1:sub_sampling/2),1)/1e3,':x','DisplayName',sprintf('Least Squares (%d iterations)',plotIndex(floor(sub_sampling/2))),'LineWidth',1,'MarkerSize',4);
set(gca,'ColorOrderIndex',3);

% Label Solutions
text(emitter_init(1)/1e3+2,emitter_init(2)/1e3,'Initial Guess','FontSize',8);
text(mean(est_plot(1,1:2,1))/1e3+5,mean(est_plot(2,1:2,1))/1e3,'Least Squares Solution','FontSize',8);
xlabel('Cross-range [km]');ylabel('Down-range [km]');
legend('Location','NorthWest');

% Plot zoomed geometry
fig_geo_b = figure;
plot(pos_emitter(1,:)/1e3,pos_emitter(2,:)/1e3,'k^','DisplayName','Transmitter');
hold on;
plot(est_plot(1,1:floor(sub_sampling/2),1)/1e3,est_plot(2,1:floor(sub_sampling/2),1)/1e3,'--x','DisplayName','LS');
plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'k','LineWidth',.5,'DisplayName','90% Error Ellipse');
xlabel('Cross-range [km]'); ylabel('Down-range [km]');
legend('Location','best');
grid off;
wd = 1.4*max(abs(crlb_ellipse(1,:) - pos_emitter(1)));
ht=wd*7/8;
ylim([-1 1]*ht/1e3 + pos_emitter(2)/1e3);
xlim([-1 1]*wd/1e3 + pos_emitter(1)/1e3);

% Plot Error
fig_err = figure;
plot(1:iters,cep50,'DisplayName','LS');
hold on;
plot(1:iters,crlb_cep50*ones(1,iters),':','DisplayName','CRLB CEP_{50}');
xlabel('Iteration Number'); ylabel('Position Error [km]');
% Set Scale
set(gca,'yscale','log');
set(gca,'xscale','log');
% Annotate Plot
text(4,6,'Least Square','FontSize',10);
text(1.5,1.2,'CRLB CEP_{50}','FontSize',10);
fprintf('The emitter is located at [%2.5f, %2.5f] km',est_iter(1,end)/1e3,est_iter(2,end)/1e3)
%% Compute and plot Isochrones 
% sensor_pos{1} = [0;0];
% for n = 2:nSensors
%     sensor_pos{n} = [xs(n);ys(n)];
%     r{n} = utils.rng(sensor_pos{n},pos_emitter);
% end
% for m = 1:nSensors
%     if abs(nSensors-m)>=2
%         xiso{m} = tdoa.drawIsochrone(sensor_pos{m},sensor_pos{m+1},r{m}-r{m+1},1000,50e3);
%     end
% end
% 
% % Plotting
% 
% figure;hold on;
% % Isochrones
% plot(xiso1(1,:),xiso1(2,:),'k:','DisplayName','Isochrone');
% hiso2=plot(xiso2(1,:),xiso2(2,:),'k:');
% hiso3=plot(xiso3(1,:),xiso3(2,:),'k:');
% utils.excludeFromLegend(hiso2);
% % Isochrone Labels
% text(mean([pos_sensor(1),pos_sensor(1)]),mean([pos_sensor(2),pos_sensor(2)])-.2,'TDOA_{1,2}');
% text(mean([pos_sensor(1),pos_sensor(1)])+.3,mean([pos_sensor(2),pos_sensor(2)]),'TDOA_{2,3}');
% plot(pos_emitter(1),pos_emitter(2),'k^','MarkerSize',8,'DisplayName','Transmitter');