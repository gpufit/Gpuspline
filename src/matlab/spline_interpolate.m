function [interpolated_data, time] = spline_interpolate(data, x, y, z)
% Interpolates 1-3D data at positions x, y, z. Internally calls
% spline_coefficients first to create a spline representation and then 
% spline_value to obtain interpolated values. If repeated interpolation
% of the same data set is needed, it's better to call spline_coefficients
% only once and then use spline_value instead of this function.

% data is single
assert(isa(data, 'single'), 'Type of data is not single');

tic;
switch nargin
    case 2
        n_dimensions = 1;
        assert(sum(size(data) > 1) == n_dimensions,'wrong input data dimensions');
        assert(isa(x, 'single'), 'Type of x is not single');
        interpolated_data = spline_interpolateMex(data, n_dimensions, x);
    case 3
        n_dimensions = 2;
        assert(sum(size(data) > 1) == n_dimensions,'wrong input data dimensions');
        assert(isa(x, 'single'), 'Type of x is not single');
        assert(isa(y, 'single'), 'Type of y is not single');
        interpolated_data = spline_interpolateMex(data, n_dimensions, x, y);
        interpolated_data = reshape(interpolated_data, numel(x), numel(y));
    case 4
        n_dimensions = 3;
        assert(sum(size(data) > 1) == n_dimensions,'wrong input data dimensions');
        assert(isa(x, 'single'), 'Type of x is not single');
        assert(isa(y, 'single'), 'Type of y is not single');
        assert(isa(z, 'single'), 'Type of z is not single');
        interpolated_data = spline_interpolateMex(data, n_dimensions, x, y, z);
        interpolated_data = reshape(interpolated_data, numel(x), numel(y), numel(z));
    otherwise
        error('wrong number of input arguments');
end
time = toc;

end