function [coefficients, time] = spline_coefficients(data)
% Calculates spline coefficients from given 1-3D data

% Wrapper around the spline mex file.

%% type checks

% data is single
assert(isa(data, 'single'), 'Type of data is not single');

%% dimensions
n_points_x = size(data, 1);
n_points_y = size(data, 2);
n_points_z = size(data, 3);

if n_points_y == 1
    n_intervals(1) = (n_points_x-1);
    n_dimensions = 1;
elseif n_points_z == 1
    n_intervals(1) = (n_points_x-1);
    n_intervals(2) = (n_points_y-1);
    n_dimensions = 2;
else
    n_intervals(1) = (n_points_x-1);
    n_intervals(2) = (n_points_y-1);
    n_intervals(3) = (n_points_z-1);
    n_dimensions = 3;
end

n_coefficients_per_interval = power(4, n_dimensions);

%% run spline taking the time
tic;
coefficients = spline_coefficientsMex(data, n_points_x, n_points_y, n_points_z, n_dimensions);
time = toc;

% reshape the output coefficients array to have dimensions
% (n_coefficients, n_intervals(1), n_intervals(2), n_intervals(3))
coefficients = reshape(coefficients, [n_coefficients_per_interval n_intervals]);

end
