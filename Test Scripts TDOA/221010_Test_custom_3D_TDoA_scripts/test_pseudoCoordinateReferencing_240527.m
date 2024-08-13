clc, clear, close all

% Define the bounds
x_bounds = [-100, 100];
y_bounds = [-100, 100];

% Sample data: create a matrix v of size m by n (for example, 200 by 200)
m = 200;
n = 200;
v = peaks(m);  % Just for example purposes, replace with your actual data
tic
% Create linearly spaced vectors for x and y
x = linspace(x_bounds(1), x_bounds(2), n);  % n points between -100 and 100 for x
y = linspace(y_bounds(1), y_bounds(2), m);  % m points between -100 and 100 for y

% Check the grid spacing
dx = (x_bounds(2) - x_bounds(1)) / (n - 1);
dy = (y_bounds(2) - y_bounds(1)) / (m - 1);

% Display the spacing for verification
fprintf('Grid spacing in x: %.4f\n', dx);
fprintf('Grid spacing in y: %.4f\n', dy);

% Coordinates where you want to find the interpolated value
x_query = 25;
y_query = 10;

% Convert coordinates to indices
x_index = round((x_query - x_bounds(1)) / dx) + 1;
y_index = round((y_query - y_bounds(1)) / dy) + 1;

% Display the indices for verification
fprintf('x_index: %d\n', x_index);
fprintf('y_index: %d\n', y_index);

% Ensure indices are within bounds
if x_index >= 1 && x_index <= n && y_index >= 1 && y_index <= m
    v_query = v(y_index, x_index);
    fprintf('The value at (x, y) = (25, 10) is: %.4f\n', v_query);
else
    error('Query point is out of bounds.');
end
toc
% Plot the original data and the indexed point for visual inspection
figure;
waterfall(x, y, v);
hold on;
plot3(x_query, y_query, v_query, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('Surface plot with direct index point');
xlabel('X-axis (km)');
ylabel('Y-axis (km)');
zlabel('v values');
hold off;
