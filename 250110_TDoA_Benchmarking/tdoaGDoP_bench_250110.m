%% Monte-Carlo Simulation for benchmarking the performance of 2D TDoA Algo 
%% against the Geometric Dilution of Precision (GDoP) (and other metrics)

% 

% Author: Abeer Nasir Chaudhry®
% 250110

clear;clc;close all

%% Initialize simulation parameters

% Initialize sensor attributes

sensorPos = [00e3, 10e3 10e3 -25e3; 
             00e3  20e3 -20e3, 00e3;
             00e3   00e3 00e3  00e3 ]; 

nSensor = size(sensorPos,2);

% Central Reference System Pairing
tdoaPairs(:,2) = (2:(nSensor))';tdoaPairs(:,1) = 1;

% Set the RMS error in time
timingErr = 30e-9; 

% Initialize source attributes
emitterPos = [-50e3;10e3;00];

nEmitter = size(emitterPos,2);

% Set the seed for the rng
rng(1,"twister");

%% Initialize the baseline hyperbolic equations as inline functions 
syms x1 x2 x3 x4 xx y1 y2 y3 y4 yy z1 z2 z3 z4 zz c tdoa12 tdoa13 tdoa14
eqn1 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x2-xx)^2 + (y2-yy)^2 + (z2-zz)^2)  == c*(tdoa12);
eqn2 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x3-xx)^2 + (y3-yy)^2 + (z3-zz)^2)  == c*(tdoa13);
eqn3 = sqrt((x1-xx)^2 + (y1-yy)^2 + (z1-zz)^2) - sqrt((x4-xx)^2 + (y4-yy)^2 + (z4-zz)^2)  == c*(tdoa14);

[x1,y1,z1]       = deal(sensorPos(1,1),sensorPos(2,1),sensorPos(3,1));
[x2,y2,z2]       = deal(sensorPos(1,2),sensorPos(2,2),sensorPos(3,2));
[x3,y3,z3]       = deal(sensorPos(1,3),sensorPos(2,3),sensorPos(3,3));
[x4,y4,z4]       = deal(sensorPos(1,4),sensorPos(2,4),sensorPos(3,4));

c = physconst('Lightspeed');
 
% Solve partial equation and create inline functions
sol = solve([eval(eqn1) eval(eqn2) eval(eqn3)],[xx yy zz]);

solHyp_x1 = matlabFunction(sol.xx(1));solHyp_x2 = matlabFunction(sol.xx(2));
solHyp_y1 = matlabFunction(sol.yy(1));solHyp_y2 = matlabFunction(sol.yy(2));
solHyp_z1 = matlabFunction(sol.zz(1));solHyp_z2 = matlabFunction(sol.zz(2));

%% Initiate MC 
MC = 10e3;

% Pre-allocate
[sol_x,sol_y,pos_est] = deal(zeros(2,MC));
dilution = zeros(1,MC);

for mc = 1:MC

    % ToA/TDoA Generator function. Generate ToA/TDoAs with set std dev and seed
    % for RNG
    [reported_tdoa, ranges] = tdoaGen_fn(sensorPos,emitterPos,timingErr,tdoaPairs);
    
    % Evaluate the referenced hyperbolic equations to get the TDoA
    % solutions
     sol_x(1:2,mc) = [solHyp_x1(reported_tdoa(1),reported_tdoa(2),reported_tdoa(3)); solHyp_x2(reported_tdoa(1),reported_tdoa(2),reported_tdoa(3))];
     sol_y(1:2,mc) = [solHyp_y1(reported_tdoa(1),reported_tdoa(2),reported_tdoa(3)); solHyp_y2(reported_tdoa(1),reported_tdoa(2),reported_tdoa(3))];
     pos_est(1:2,mc) = [sol_x(2,mc).' sol_y(2,mc).'];

    % Compute dilution at the emitter position 
    [~,dilution(mc)] = modfn_gdop2D(sensorPos,emitterPos,timingErr);

end

% Compute the position error covariance matrix using the mean of estimated
% positions. Mean is used because it corresponds more directly to th
% precision of the system and is used for benchmarking against the CRLB
meanPosEst = mean(pos_est,2);
posErrCovMtx = ((pos_est - meanPosEst)*(pos_est - meanPosEst).')/MC;

% Get the covariance matrix for the truePosition of the emitter to create
% the error ellipse for benchmarking
[FIM,~] = modfn_gdop2D(sensorPos,emitterPos,timingErr);

%% Performance Evaluation
range_err = vecnorm(pos_est - emitterPos(1:2,:)); % Compute the error in position
meanRangerr= mean(range_err);                     % Compute the median of the position error
ellipseMed = median(pos_est,2);                   % Identify the median of the error ellipse
ellipseMedErr = vecnorm(ellipseMed - emitterPos(1:2,:)); % Compute the positional error of the median of the dist vs true pos

% Bearing Error
bearTrue = atan2d(emitterPos(2,1) - (sensorPos(2,1)),(emitterPos(1,1) - sensorPos(1,1))); % True bearing w.r.t CRS
bearEst = atan2d((pos_est(2,:) - sensorPos(2,1)),(pos_est(1,:) - sensorPos(1,1))); % Bearing of ests w.r.t CRS
bearErr = bearTrue - bearEst;
meanbearErr = mean(bearErr);

% Error Ellipses
crlb_ellipse = drawErrorEllipse(emitterPos(1:2,:),inv(FIM),100,90); % Create the error ellipse representing the best theoretical performance of an unbiased estimator
system_ellipse = drawErrorEllipse(emitterPos(1:2,:),posErrCovMtx,100,90); % Create the error ellipse for the current system performance

% CEP50
cep = CEP50(posErrCovMtx);                       % Compute the CEP50 for the system based on pos err cov mtx
cepCRLB = CEP50(inv(FIM));                       % Compute the theoretical lower bound on the CEP50

% Check what percentage of TDoA solutions comply with the CRLB via the
% CEP50
cepcrlbCheck = length(range_err(range_err<cepCRLB))*100/MC;
cepsystemCheck = length(range_err(range_err<cep))*100/MC;

fprintf(['Mean Positional Error: %0.3f m\n' ...
    'Positional Error from Median of Ellipse: %0.3f m\n' ...
    'System CEP50: %0.3f m\n' ...
    'Theoretical Min CEP50: %0.3f m\n' ...
    'CEP50 compliance Percentage (CRLB): %0.3f %%\n', ...
    'CEP50 compliance Percentage (System): %0.3f %%\n'], ...
    meanRangerr, ellipseMedErr, cep,cepCRLB,cepcrlbCheck,cepsystemCheck)

%%  Visualize results
figure; hold on; grid on;
sen = scatter(sensorPos(1,:)/1e3,sensorPos(2,:)/1e3,'filled','o','DisplayName','Sensor Network','MarkerEdgeColor','k','MarkerFaceColor','b');
est = scatter(pos_est(1,:)/1e3,pos_est(2,:)/1e3,'.','DisplayName','TDoA Solutions','MarkerEdgeColor','k','MarkerFaceColor','y');
med = scatter(ellipseMed(1,:)/1e3, ellipseMed(2,:)/1e3,'+','DisplayName','Median of the Distribution','MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',2);
em = scatter(emitterPos(1,:)/1e3,emitterPos(2,:)/1e3,'filled','^','DisplayName','True Emitter Position','MarkerEdgeColor','k','MarkerFaceColor','r');
plot(crlb_ellipse(1,:)/1e3,crlb_ellipse(2,:)/1e3,'g','LineWidth',.5,'DisplayName','90% Error Ellipse (CRLB)');
plot(system_ellipse(1,:)/1e3,system_ellipse(2,:)/1e3,'color',[1 0.7529 0],'LineWidth',.5,'DisplayName','90% Error Ellipse (System)');

%% Compute the eigen vectors 
covMtx = cov(pos_est(1,:)/1e3,pos_est(2,:)/1e3);

[eigVec ,eigVal] = eig(covMtx);
d = sqrt(diag(eigVal));
hold on
eig1plt1 = quiver(emitterPos(1)/1e3,emitterPos(2)/1e3,eigVec(1,1),eigVec(2,1),d(1),'Color','g','DisplayName','Eigen Vector','HandleVisibility','off');
eig1plt2 = quiver(emitterPos(1)/1e3,emitterPos(2)/1e3,eigVec(1,2),eigVec(2,2),d(2),'Color','g','DisplayName','Eigen Vector');

xlabel('Cross-Range (km)');ylabel('Down-Range (km)')
legend('Location','northeastoutside')

figure;
plot(1:length(range_err),range_err/1e3,'DisplayName','Positional error in TDoA measurements');
hold on
plot(1:length(range_err),dilution/1e3,'DisplayName','Dilution of Precision')
plot(1:length(range_err),meanRangerr/1e3*ones(1,length(range_err)),'DisplayName','Mean of TDoA Measurement Errors','LineWidth',3)
plot(1:length(dilution),mean(dilution)/1e3*ones(1,length(dilution)),'DisplayName','Mean of DoP','LineWidth',3)
xlabel('Monte-Carlo Runs');ylabel('Displacement from True Emitter Position (km)')
legend;

figure;
plot(1:length(bearErr), bearErr,'DisplayName','Error in bearing'); hold on
plot(1:length(bearErr),(meanbearErr)*ones(1,length(bearErr)),'DisplayName','Mean of Bearing Error')
xlabel('Monte-Carlo Runs');ylabel('Error in Bearing (deg')
legend;

% 
% % LOG
% 1. Date created (updated): 250110 (250114)

% 2. Added the error ellipse calculated using the CRLB computed using the
% modfn_gdop2D.

% 3. The CEP50 is a scalar measure of uncertainty, representing the radius 
% of a circle that contains 50% of the position estimates. It simplifies 
% the 2D error distribution into a single value. CEP50 computed using the
% positional error covariance matrix represents the uncertainty derived
% from the observed spread o the position estimates (in a Monte Carlo Sim).
% It is the actual peroformance of the estimator including all biases and
% inefficiences.
%  The CEP50 calculated using the CRLB is the theoretical minimum error
%  that can be achieved by an unbiased estimator for the given scenario. It
%  can be used as a benchmark to  evaluate the estimators performance.

% 4. The error ellipse represents the 2D uncertainty in the position 
% estimates based on the error covariance matrix. It visually shows the 
% confidence region for the position estimates.The semi-major and 
% semi-minor axes of the ellipse are derived from the eigenvalues of 
% the covariance matrix, which quantify the variances in the principal 
%  directions of the error distribution.

% Both CEP50 and error ellipse are interconnected and complementary, not 
% independent, in evaluating estimator precision.

% 5. Bearing Error calculated and added to the performance evaluation

% 6. It was considered to save the state of the RNG and use that to supply
% the same randomization to the ToA generation and to the calculation of
% the CRLB/GDoP. When this was done, it was observed that the theoretical 
% minumum CEP50 derived using the CRLB and the CRLB error ellipse both
% the statistical metric were smaller than those exhibited by the system,
% which is unrealistic incoherent. After careful literatur review, it was 
% determined that to ensure the accuracy and validity of the performance
% evaluation, it is crucial to maintain independence between the
% randomization of the ToAs and of the covariance matrix used to build the
% CRLB. Reusing the same random number generator (RNG) state for both 
% processes introduces unintended correlations between the measurement 
% noise and the error covariance matrix for the GDoP/CRLB, leading to 
% misleadingly improved performance metrics. The system error can never be 
% lower than the theoretical minimum error computed using the CRLB. The 
% noise in the TDoAs is now partially accounted for by the variance in Q.
% This correlation reduces the uncertainty in the localization solution, 
% leading to an artificially improved performance (smaller CEP50 and 
% tighter error ellipses).

