function [X,Y,GDoP] = gdop2Dsd(sensorPos,timingErr,res,gridSpan,refSen)

%% 2D GDoP calculation

% Range difference based centralized reference pairing calculation of GDoP.
% The Error Cov Mtx that is used is based on the following paper: 
% Lowest GDOP in 2-D scenarios N.Levanon

% INPUTS: 
        % 1. sensorPos: 2xN matrix of N sensor positions in 2D (matrix)
        % 2. res: Resolution in meters for manifestation of candidate 
        % emitter position (scalar)
        % 3. timingErr: RMS value of timing measurement error (scalar)
        % 4. gridSpan: A [2 x 2] matrix containing the upper and lower
        % bounds of the grid axes. Columns contain the upper and lower
        % bound for x-axis, while rows hold the bounds for y-axis. (matrix)
        % 5. refSen: Reference sensor index (scalar)

% OUTPUT(S):
        % 1. X: 2D grid coordinates for x-axis.  (matrix)
        % 2. Y: 2D grid coordinates for y-axis.  (matrix)
        % 3. GDoP: Geometric Dilution of Precision (matrix)
        

%% Process
% 1. Initialize parameters
% 2. Create candidate emitter positions
% 3. Create the Jacobian for each position by computing the range and
% directional derivatives
% 4. Update the Fischer Information Matrix (FIM)
% 5. Compute and store the GDoP for each position by taking the sqrt of 
% trace of the FIM

% Abeer Chaudhry©
% Jan 2025

%% Checks
if size(sensorPos,1) ~= 2
    error('Input sensor positions do not match the required matrix dimensions (2 x Nsensors)');
end

if nargin < 5  || isempty(refSen)
    refSen = 1;
end

if nargin < 4 
    gridSpan = [-400e3,400e3;-400e3, 400e3];
end

if nargin < 3 
    res = 2e3;
end

%% Initialize Parameters
nSensors = size(sensorPos,2);
sensVec = setdiff(1:nSensors,refSen); % Parse vector for other sensors

% Restructure the covariance matrix based on the choice of reference sensor
ErrCov = fn_RDoaErrCov(nSensors,refSen,timingErr);

InvErrCov = pinv(ErrCov); % Compute the psuedo-inverse of the rectangular, non-singular matrix.
% The structure of the covariance matrix is determined by the the
% case/mechanism that is employed for calculating the Jacobian for the CRLB
% which determines the GDoP of a network. Here the non-diagonal elements
% are 1 as opposed to being 0 which represents the statistical dependencies
% between the sensors with the reference sensor as a common delay 
% measurement used to generate the differences. The 2 in the diagonal 
% reflect that the range-difference measurements are a relative delay. More
% in Log.

% Create Candidate Emitter Position Grid
Xspan = gridSpan(1,:);
Yspan = gridSpan(2,:);

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
        
        for n = 1:length(sensVec)
            H(n,:) = [(emitterPos(1) - sensorPos(1,sensVec(n)))/range(sensVec(n)) - (emitterPos(1) - sensorPos(1,refSen))/range(refSen);
                        (emitterPos(2) - sensorPos(2,sensVec(n)))/range(sensVec(n)) - (emitterPos(2) - sensorPos(2,refSen))/range(refSen)];
        end
        
        FIM = inv(H.'*InvErrCov*H);
        GDoP(i, j) = sqrt(trace(FIM));
    end
end

%% Log:
% 1. In range-difference measurements, we measure the difference in ranges 
% between pairs of sensors rather than the absolute ranges. The measurement
% noise in one sensor contributes twice in the difference (once positively 
% and once negatively in the two range terms).
% This creates correlated errors between different range-difference 
% measurements.

% var(delta_r) = var(ri) + var(rj) = sig^2 + sig^2 =2*sig^2
% Additionally, two different range differences (e.g., 𝑟_i ), leading to
% correlation between the two measurements which is sig^2.

% 2. There is a difference between ndgrid and meshgrid. Just FYI

% 3. A function namely 'restructCovMtx' is added to the script. Since the
% GDoP is dependent upon the sensor geometry and is calculated with
% reference to a central sensor, the choice of the reference sensor affects
% the geomtric relationship of the network to that of the emitter and thus
% changes the DoP. The script has been updated to allow the user to select
% a reference sensor index (set at 1by default). As a result of a reference
% sensor other then i=1, the covariance matrix must also be resampled with
% the first row of the (NSensor x NSensor) matrix representing the
% reference sensor.
end