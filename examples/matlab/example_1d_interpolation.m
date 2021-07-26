function example_1d_interpolation()
%%% An example for interpolating data points calling                 %%%
%%% a cubic spline interpolation routine implemented in C.           %%%
%%% The data is upsampled, cut, stretched and shifted.               %%%

% input data
y = single([0,0,0.2,1,1.1,1.3,2,2.5,3,4,4.25,4,3,2.5,2,1.3,1.1,1,0.2,0,0])';
x = single(0:numel(y)-1);
center_index = x(end)/2;

% interpolation parameters
edge = 1.4;
width = 1.1;
shift = 1.2;
sampling_factor = 0.5;

% interpolation
xq = x(1):sampling_factor:x(end);
xq = xq(xq >= edge & xq <= max(xq)-edge);
xq = xq / width;
xq = xq+(center_index*(1-1/width))-shift;
yq = spline_interpolate(y,xq);

% figure
plot(x,y,'-bs');
hold on;
xx = linspace(x(1)+edge,x(end)-edge,numel(xq));
plot(xx,yq,'-rx');
xlim([0 20]);
ylim([0 max(y)+max(y)/10]);
legend('original', 'interpolated');
hold off;

end