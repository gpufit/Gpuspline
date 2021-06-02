function [values, time] = spline_values(coefficients, x, y, z)
% calculates spline values given a spline representation and x, y and z
% values
%
% y and z parameter is optional (only if you want 2D/3D values)
%
% note: x, y, z are a range of values along their coordinate axis and the total
% number is the cartesian product of them all!
%
% TODO this should change, we would rather give full x, y, z arrays
% allowing arbitrary combinations of x, y, z values because internally
% anyway the intervalls are computed in which each x,y,z value falls

%% type checks
assert(isa(coefficients, 'single'),'Type of spline coefficients is not single');
assert(isa(x, 'single'),'Type of x is not single');
if exist('y','var'); assert(isa(y, 'single'),'Type of y is not single'); end
if exist('z','var'); assert(isa(z, 'single'),'Type of z is not single'); end

%% dimensions check
switch nargin
    case 2
        n_dimensions = 1;
        y = single(1);
        z = single(1);
    case 3
        n_dimensions = 2;
        z = single(1);
    case 4
        n_dimensions = 3;
    otherwise
        error('wrong number of input arguments')
end

dimensions = zeros(5,1);
for i=1:5; dimensions(i) = size(coefficients,i); end

if dimensions(5) > 1
    assert(sum(dimensions > 1) == n_dimensions + 2,...
        'Wrong spline coefficients dimensions.');
else
    assert(sum(dimensions > 1) == n_dimensions + 1,...
        'Wrong spline coefficients dimensions');
end

tic;

coefficients = reshape(coefficients, [], dimensions(5));

values = zeros(numel(x)*numel(y)*numel(z),dimensions(5),'single');

for ch = 1:dimensions(5)
    if n_dimensions == 1
        values(:,ch) = spline_valuesMex(...
            coefficients(:,ch),...
            n_dimensions,...
            dimensions(2),...
            dimensions(3),...
            dimensions(4),...
            x);
    elseif n_dimensions == 2
        values(:,ch) = spline_valuesMex(...
            coefficients(:,ch),...
            n_dimensions,...
            dimensions(2),...
            dimensions(3),...
            dimensions(4),...
            x, y);
    elseif n_dimensions == 3
        values(:,ch) = spline_valuesMex(...
            coefficients(:,ch),...
            n_dimensions,...
            dimensions(2),...
            dimensions(3),...
            dimensions(4),...
            x, y, z);
    end
end

time = toc;

values = reshape(values,numel(x),numel(y),numel(z),dimensions(5));

end
