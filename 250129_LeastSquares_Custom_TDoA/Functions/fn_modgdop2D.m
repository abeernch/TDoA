function [FIM,GDoP] = fn_modgdop2D(sensorPos,emitterPos,timingErr,refSen)

%% 2D GDoP calculation

% Range difference based centralized reference pairing calculation of GDoP.
% The Error Cov Mtx that is used is based on the following paper: 
% Lowest GDOP in 2-D scenarios N.Levanon

% This is a modified version of the orignal function where instead of a
% grid, a single emitter position is provided to compute the GDoP at that
% position for a range of different errors in measurement. 
% The purpose of this modified function is to benchmark the performance of
% the TDoA algorithm by performing Monte-Carlo simulations.

% INPUTS: 
        % 1. sensorPos: 2xN matrix of N sensor positions in 2D (matrix)
        % 2. emitterPos: 2xN matrix of emitter positions in 2D (matrix) 
        % 3. timingErr: RMS value of timing measurement error (scalar)
        % 4. refSen: Reference sensor index (scalar)

% OUTPUT(S):
        % 1. FIM: Fischer Information Matrix (matrix)
        % 2. GDoP: Geometric Dilution of Precision (matrix)
        

%% Process
% 1. Initialize parameters
% 2. Create the Jacobian for each position by computing the range and
% directional derivatives
% 4. Update the Fischer Information Matrix (FIM)
% 5. Compute and store the GDoP for the emitter position by taking the sqrt of 
% trace of the FIM

% Abeer Chaudhry©
% Jan 2025

%% Checks
sensorPos = sensorPos(1:2,:);
emitterPos = emitterPos(1:2,:);
if size(sensorPos,1) ~= 2 || size(emitterPos,1) ~= 2
    error('Input positions do not match the required matrix dimensions (2 x Nsensors)');
end

if nargin < 4  || isempty(refSen)
    refSen = 1;
end

%% Initialize Parameters
nSensors = size(sensorPos,2);
sensVec = setdiff(1:nSensors,refSen); % Parse vector for other sensors

% Restructure the covariance matrix based on the choice of reference sensor
ErrCov = fn_RDoaErrCov(nSensors,refSen,timingErr);

InvErrCov = inv(ErrCov); % Compute the psuedo-inverse of the rectangular, non-singular matrix.
% The structure of the covariance matrix is determined by the the
% case/mechanism that is employed for calculating the Jacobian for the CRLB
% which determines the GDoP of a network. Here the non-diagonal elements
% are 1 as opposed to being 0 which represents the statistical dependencies
% between the sensors with the reference sensor as a common delay 
% measurement used to generate the differences. The 2 in the diagonal 
% reflect that the range-difference measurements are a relative delay. More
% in Log.

warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning

%% Repeat CRLB for each of the candidate emitter positions

range = vecnorm(emitterPos - sensorPos);

% Initialize the jacobian matrix
H = zeros(nSensors - 1, 2);

for n = 1:length(sensVec)
    H(n,:) = [(emitterPos(1) - sensorPos(1,sensVec(n)))/range(sensVec(n)) - (emitterPos(1) - sensorPos(1,refSen))/range(refSen);
        (emitterPos(2) - sensorPos(2,sensVec(n)))/range(sensVec(n)) - (emitterPos(2) - sensorPos(2,refSen))/range(refSen)];
end

FIM = (H.'*InvErrCov*H);      % FIM
GDoP = sqrt(trace(inv(FIM))); % Sqrt of the trace of CRLB (inv of FIM)

%% Log: 
% date created (updated): 250107 (250114)
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

% 4. Another output has been added to the function. It now also outputs the
% FIM. This is sometimes required to get the CRLB for a position.
% 5. Inverse used instead of pseudo inverse since the matrices are expected
% to be well-conditioned.

end