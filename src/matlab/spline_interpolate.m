function [interpolated_data, time] = spline_interpolate(data, x, y, z)
% interpolates data at positions x, y, z by using a combination of
% spline_coefficients and spline_value

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