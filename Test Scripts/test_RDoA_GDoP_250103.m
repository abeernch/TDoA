%% 2D GDoP calculation

% Range difference based centralized reference pairing calculation of GDoP.
% The Error Cov Mtx that is used is based on the following paper: 
% Lowest GDOP in 2-D scenarios N.Levanon


%% Process
% 1. Initialize parameters
% 2. Create candidate emitter positions
% 3. Create the Jacobian for each position by computing the range and
% directional derivatives
% 4. Update the Fischer Information Matrix (FIM)
% 5. Compute and store the GDoP for each position by taking the sqrt of 
% trace of the FIM

% Abeer Chaudhry ©
% Jan 2025
clc;clear;
%% Initialize Parameters
c = physconst('LightSpeed');

sensorPos = [00e3	-20e3	10e3	10e3
             0	  0 	   20e3	 -20e3];

nSensors = size(sensorPos,2);

timingErr = 1e-9;

InvErrCov = pinv(timingErr^2*(eye(nSensors-1)+1).*c^2); % Compute the psuedo-inverse of the rectangular, non-singular matrix.
% The structure of the covariance matrix is determined by the the
% case/mechanism that is employed for calculating the Jacobian for the CRLB
% which determines the GDoP of a network. Here the non-diagonal elements
% are 1 as opposed to being 0 which represents the statistical dependencies
% between the sensors with the reference sensor as a common delay 
% measurement used to generate the differences. The 2 in the diagonal 
% reflect that the range-difference measurements are a relative delay. More
% in Log.

% Create Candidate Emitter Position Grid
Xspan = [-400e3,400e3];
Yspan = [-400e3,400e3];
res = 2e3;              % Grid resolution in meters

% A grid of potential emitter locations to be evaluated is created. The
% bounds are defined by XSpan and YSpan with a resolution of 'res' which is
% translated into the equivalent no. of points for the grid using:
% pts = (XSpan(2) - XSpan(1))/res

[X,Y] = meshgrid(linspace(Xspan(1),Xspan(2),(Xspan(2) - Xspan(1))/res),linspace(Yspan(1),Yspan(2),(Yspan(2) - Yspan(1))/res));

GDoP = zeros(size(X)); % Fischer Info Matrix (Dilution matri) preallocation

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning

%% Repeat CRLB for each of the candidate emitter positions
for i = 1:size(X,1)
    for j = 1:size(X,2)
        emitterPos = [X(i,j),Y(i,j)];
        range = vecnorm(emitterPos.' - sensorPos);
        
        % Initialize the jacobian matrix
        H = zeros(nSensors - 1, 2);
        
        for n = 2:nSensors
            H(n-1,:) = [(emitterPos(1) - sensorPos(1,n))/range(n) - (emitterPos(1) - sensorPos(1,1))/range(1);
                        (emitterPos(2) - sensorPos(2,n))/range(n) - (emitterPos(2) - sensorPos(2,1))/range(1)];
            
        end
        
        FIM = inv(H.'*InvErrCov*H);
        GDoP(i, j) = sqrt(trace(FIM));
    end
end

contourLevelsxy = 100:100:2000;
figure;
scatter(sensorPos(1,:),sensorPos(2,:))
hold on
contour(X,Y,GDoP,contourLevelsxy,'LineStyle','--','LineColor','k')



%% Log:
% 1. In range-difference measurements, we measure the difference in ranges between pairs of sensors rather than the absolute ranges.
% The measurement noise in one sensor contributes twice in the difference (once positively and once negatively in the two range terms).
% This creates correlated errors between different range-difference measurements.
% var(delta_r) = var(ri) + var(rj) = sig^2 + sig^2 =2*sig^2
% Additionally, two different range differences (e.g., 𝑟_i ), leading to
% correlation between the two measurements which is sig^2.
% 2. There is a difference between ndgrid and meshgrid. Just FYI
