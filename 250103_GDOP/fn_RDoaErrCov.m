function ErrCov = fn_RDoaErrCov(nSensors, refSen,timingErr)

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
            % 1. nSensors: no of sensors in geometry (scalar)
            % 2. refSen: Reference sensor index (scalar)
            % 3. timinErr: RMS value of timing measurement error (scalar)
    
    % OUTPUTS:
            % 1. ErrCov: RDoA Error Covariance Matrix (nSensor-1 x nSensor-1)

    % Abeer Chaudhry ©
    % January 2025

    c = physconst('LightSpeed');
    
    % Original Covariance Matrix
    Q = eye(nSensors)*timingErr^2*c^2;
    
    % Create Difference Operator Matrix
    sensor_pairs = setdiff(1:nSensors, refSen);
    D = sparse(length(sensor_pairs), nSensors);
    
    for i = 1:length(sensor_pairs)
        D(i, sensor_pairs(i)) = 1;
        D(i, refSen) = -1;
    end
    
    % Reduce the Covariance Matrix
    ErrCov = D*Q*D.';
    

end