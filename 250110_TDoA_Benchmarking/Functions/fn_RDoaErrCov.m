function ErrCov = fn_RDoaErrCov(nSensor, refSen,timingErr)

    % Since the calculations are based on the RDoA/TDoA relative to a
    % reference sensor, the matrix and dependence in the variances vary
    % accordingly. For example, if i=3 was selected as a reference sensor,
    % the covariance matrix must be resampled accordingly with the order of
    % the rows updated accordingly and the difference now taken w.r.t the
    % 3rd sensor as the reference sensor.

    % This is carried out using the Difference Operator Matrix which
    % performs the range difference conversion in vectorized form relative
    % to the selected reference sensor. The general covariance matrix Q
    % represents the noise noise correlation between sensor ranges which
    % needs to be converted to RDoA based noise correlation. This is
    % carried out using the following linear transformation:

    % Q_RDoA = D*R*D.'

    % D*R: Transforms the covariance matrix into the difference space.
    % D.' : Maps it back to ensure the resulting matrix represents the 
    % correlation between range differences.

    % This process effectively translates the following iterative 
    % transformation:
    
    %    [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
    %       where:  ai, aj are the i-th and j-th reference indices
    %               bi, bj are the i-th and j-th test indices
    
    % INPUTS: 
            % 1. nSensor: no of sensors in geometry (scalar)
            % 2. refSen: Reference sensor index (scalar)
            % 3. timinErr: RMS value of timing measurement error (scalar)
    
    % OUTPUTS:
            % 1. ErrCov: RDoA Error Covariance Matrix (nSensor-1 x nSensor-1)

    % Abeer Chaudhry ©
    % January 2025

    c = physconst('LightSpeed');
    
    % Original Covariance Matrix
    % Q = eye(nSensor)*timingErr^2*c^2;
    
    % To randomize the covariance for the MC simulation using a Gaussian
    % dist rand number, it is important to note that the positive
    % semi-definite-ness of Q must be maintained. If randn^2 is multiplied
    % with Q to add randomness, it must be noted that Gaussian distribution
    % when squared, becomes a Chi-squared distribution and can skew the
    % randomness.
    % The correct way would be to use:
    % Q = R*R' where R = randn(Nsensor)
    % Here, since this is the cov mat for the abs range measurements, the
    % off-diagonal terms are all 0, indicating no correlation (as it should
    % be). Therefore, to maintain the structure followed by its
    % restructuring with the RDoA cov mat linear tfm, the diagonal terms
    % are randomized. They are randomized non-uniformly/uniformly to the 
    % diagonal terms to represent independent variances.
    
    % Non-uniformly randomized variances for heterogenous sensor network
    % Q = eye(nSensor).*abs(randn(1,4))*timingErr^2*c^2;

    % Uniformly randomized variances for homogenous sensor network
    Q = eye(nSensor).*abs(randn(1))*timingErr^2*c^2;
    
    % Create Difference Operator Matrix
    sensor_pairs = setdiff(1:nSensor, refSen);
    D = zeros(length(sensor_pairs), nSensor);
    
    for i = 1:length(sensor_pairs)
        D(i, sensor_pairs(i)) = 1;
        D(i, refSen) = -1;
    end
    
    % Reduce the Covariance Matrix
    ErrCov = D*Q*D.';

end

% LOG: 
% 1. Created (Updated) 250111 (250115)

% 2. Choice between uniformly randomized and non-uniformly randomized 
% variances for the Err Cov matrix Q. It is observed that with uniformly
% randomized variances for Q, the CRLB and hence the Error ellipse derived
% from it to visualize the theoretical minimum bound is more consistent
% with expectations based on the specific scenario dynamics. It is
% suggested that the randomized variances are required in the case when the
% sensors are vastly heterogenous in many aspects ,timing errors, clock,
% etc and other hardware and environmental conditions. As far as the effect
% of sensor geometry on the localization performance is concerened, the
% uniform randomization suits best.

% 3. 