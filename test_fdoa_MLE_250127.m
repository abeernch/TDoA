% FDoA-Based Emitter Localization with Stationary Sensors and Moving Emitter
% This script estimates the position of a moving emitter using FDoA measurements.

% Clear workspace and close figures
clear;
clc;
close all;

% Parameters
c = 3e8; % Speed of light (m/s)
f0 = 1e9; % Emitter's carrier frequency (Hz)
sigma_m = 100; % Standard deviation of FDoA measurement noise (Hz)
area_size = 400e3; % Size of the search area (m)
grid_resolution = 500; % Grid resolution (m)
num_samples = 1e1; % Number of samples for sample-based method

% Sensor positions (stationary)
sensor_positions = [0, 0; 0, 20e3; 20e3*cosd(-30), 20e3*sind(-30); 20e3*cosd(210), 20e3*sind(210)]; % Sensor positions (m)

% True emitter position and velocity (for simulation)
true_emitter_position = [160.254e3, 39.40e3]; % True emitter position (m)
true_emitter_velocity = [100, 50]; % True emitter velocity (m/s)

% Generate FDoA measurements
% Compute true range rate differences
num_sensors = size(sensor_positions, 1);
true_range_rate_differences = zeros(num_sensors - 1, 1);

for i = 2:num_sensors
    % Range rate difference between sensor 1 and sensor i
    range_rate1 = dot(true_emitter_velocity, true_emitter_position - sensor_positions(1, :)) / norm(true_emitter_position - sensor_positions(1, :));
    range_rate2 = dot(true_emitter_velocity, true_emitter_position - sensor_positions(i, :)) / norm(true_emitter_position - sensor_positions(i, :));
    true_range_rate_differences(i - 1) = (f0 / c) * (range_rate1 - range_rate2);
end

% Add noise to the measurements
noisy_range_rate_differences = true_range_rate_differences + sigma_m * randn(size(true_range_rate_differences));

% Grid-based estimation
% Generate grid
x_grid = -area_size/2:grid_resolution:area_size/2;
y_grid = -area_size/2:grid_resolution:area_size/2;
[X, Y] = meshgrid(x_grid, y_grid);
grid_points = [X(:), Y(:)];

% Compute log-likelihood for each grid point
log_likelihood_values = zeros(size(grid_points, 1), 1);
for i = 1:size(grid_points, 1)
    % Predicted range rate differences for the current grid point
    predicted_range_rate_differences = zeros(num_sensors - 1, 1);
    for j = 2:num_sensors
        range_rate1 = dot(true_emitter_velocity, grid_points(i, :) - sensor_positions(1, :)) / norm(grid_points(i, :) - sensor_positions(1, :));
        range_rate2 = dot(true_emitter_velocity, grid_points(i, :) - sensor_positions(j, :)) / norm(grid_points(i, :) - sensor_positions(j, :));
        predicted_range_rate_differences(j - 1) = (f0 / c) * (range_rate1 - range_rate2);
    end

    % Log-likelihood (assuming Gaussian noise)
    residuals = noisy_range_rate_differences - predicted_range_rate_differences;
    log_likelihood_values(i) = -0.5 * sum(residuals.^2); % Proportional to log-likelihood
end

% Find the grid point with the maximum log-likelihood
[~, max_idx] = max(log_likelihood_values);
grid_estimates = grid_points(max_idx, :);

% Sample-based estimation
% Generate random samples
samples = area_size * (rand(num_samples, 2) - 0.5);

% Compute log-likelihood for each sample
log_likelihood_values_samples = zeros(num_samples, 1);
for i = 1:num_samples
    % Predicted range rate differences for the current sample
    predicted_range_rate_differences = zeros(num_sensors - 1, 1);
    for j = 2:num_sensors
        range_rate1 = dot(true_emitter_velocity, samples(i, :) - sensor_positions(1, :)) / norm(samples(i, :) - sensor_positions(1, :));
        range_rate2 = dot(true_emitter_velocity, samples(i, :) - sensor_positions(j, :)) / norm(samples(i, :) - sensor_positions(j, :));
        predicted_range_rate_differences(j - 1) = (f0 / c) * (range_rate1 - range_rate2);
    end

    % Log-likelihood (assuming Gaussian noise)
    residuals = noisy_range_rate_differences - predicted_range_rate_differences;
    log_likelihood_values_samples(i) = -0.5 * sum(residuals.^2); % Proportional to log-likelihood
end

% Find the sample with the maximum log-likelihood
[~, max_idx_sample] = max(log_likelihood_values_samples);
sample_estimates = samples(max_idx_sample, :);

% Display results
fprintf('True Emitter Position: (%.2f, %.2f) km\n', true_emitter_position / 1e3);
fprintf('Grid-Based Estimate: (%.2f, %.2f) km\n', grid_estimates / 1e3);
fprintf('Sample-Based Estimate: (%.2f, %.2f) km\n', sample_estimates / 1e3);