function filter = InitHighSpeedKF(detection)
% This function initializes a constant velocity Kalman filter and sets a
% higher initial state covariance on velocity components to account for
% high speed of the object vehicles. 
filter = initcvkf(detection);
filter.StateCovariance(2:2:end,2:2:end) = 500*eye(3);
end
