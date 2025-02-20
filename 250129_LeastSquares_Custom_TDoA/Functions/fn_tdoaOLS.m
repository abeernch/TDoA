function [error, posPredict] =  fn_tdoaOLS(sensorPos, refSen, pos_est, measRange, timingErr, tol, iter)
    %% Ordinary Least Squares for TDoA Solution
    % The function iteratively minimizes the mean square error between the
    % measurement and the position estimate, given and initial position
    % estimate.
    % The initial guess for the source position provided to the algorithm
    % is the one obtained using the analytical solution. 
    % The error covariance for a given Gaussian distributed noise with zero
    % mean and timingErr^2 variance is computed using the function 
    % "fn_RDoaErrCov".
    % 
    % The residual is computed as the difference between the measured RDoAs
    % and the RDoAs based on the current estimate of the LSE algo. 
    % 
    % The Jacobian is also comuted for the given estimates and the sensor
    % geometry.

    % INPUTS: 
            % 1. sensorPos: Sensor positions (m) (Ndim x Nsensor)
            % 2. refSen: Index of the reference sensor (scalar)
            % 3. pos_est: Initial estimate of source position (NDim x 1)
            % 4. measRange: The measured ranges from sensors to the source
            % (Nsensor x 1)
            % 5. timingErr: RMS timing supplied to the Err Cov function
            % 6. tol: Desired position error tolerance (stopping condition)
            % 7. iter: No of LSE iterations to perform

    % OUTPUTS:
            % 1. error: LSE error (1 x iter)
            % 2. posPredict: Iteration-by-iteration estimated source
            % position

    nSensor = size(sensorPos,2);

    % Compute the error covariance
    ErrCov = fn_RDoaErrCov(nSensor, refSen,timingErr);
    
    % Initialize the estimator parameters
    error = Inf*ones(1,iter);                     % Initialize error
    J = zeros(nSensor - 1, 2);                    % Initialize Jacobian
    [posPredict,posoffset] = deal(zeros(2,iter)); % Preallocate arrays
    sensVec = setdiff(1:nSensor,refSen);          % Sensor selector vector
    
    % Compute the RDoA from the measurements with the reference sensors
    measuredRDoA = (measRange(refSen) - measRange(sensVec));
    % measuredRDoA = measRange;
    % [xx,~] = deal([1e3;1e3]);
    
    for i = 1:iter
        % Compute the range from each sensor to the hyperbolic position estimate
        estRange = vecnorm(pos_est - sensorPos(1:2,:));
    
        % Compute the Jacobian for each update
        for n = 1:length(sensVec)
            J(n,:) = [(pos_est(1) - sensorPos(1,sensVec(n)))/estRange(sensVec(n)) - (pos_est(1) - sensorPos(1,refSen))/estRange(refSen);
                (pos_est(2) - sensorPos(2,sensVec(n)))/estRange(sensVec(n)) - (pos_est(2) - sensorPos(2,refSen))/estRange(refSen)];
        end
    
        % Compute the RDoA with the reference sensor for the current estimate
        estimateRDoA = (estRange(refSen) - estRange(sensVec)).';
    
        % Compute the residual (measuredRDoA - estimateRDoA)
        y = measuredRDoA - estimateRDoA;
        ErrCov = cov(y)*eye(3);
        % Compute the position offset
        posoffset(:,i) = (J.'*inv(ErrCov)*J+1e-6*eye(2))\(J.'*inv(ErrCov)*y);
    
        % Update prediction
        posPredict(:,i) = pos_est - posoffset(:,i);
      
        % Compile and update error
        error(:,i) = norm(posoffset(:,i));
    
        % Update variables
        pos_est = posPredict(:,i);
        prevErr = error(:,i);

        % Check if tolerance for position error reached (stopping condition)
        if error(:,i) < tol
            fprintf('Tolerance reached at iteration # %1d, breaking operation\n',i)
            break
        end
        posPredict = posPredict(:,1:i);
        error = error(1:i);
    end
end

%% LOG
% 1. Date created (updated): 250110 ()

% 2. The error term here represents the correction step size (or position 
% update magnitude) at each iteration. It is not the absolute localization 
% error but rather an indication of how much the position estimate is being
% refined in each step.