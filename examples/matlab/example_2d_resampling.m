function example_2d_resampling()

%% psf size
size_x = 15;
size_y = 20;

%% derived values
x = single(0 : size_x - 1)';
y = single(0 : size_y - 1);

x_up = single(0 : 0.1 : size_x - 1)';
y_up = single(0 : 0.1 : size_y - 1)';

x_down = single(0 : 2 : size_x - 1)';
y_down = single(0 : 2 : size_y - 1)';

%% PSF parameters
psf_parameters = single([100, (size_x-1)/2, (size_y-1)/2, 3, 10]);

%% calculate PSF
psf = calculate_psf(x, y, psf_parameters);

%% calculate spline coefficients
coefficients = spline_coefficients(psf);

%% generate upsampled psf
psf_up = spline_values(coefficients, x_up, y_up);

%% generate downsampled psf
psf_down = spline_values(coefficients, x_down, y_down);

%% figure
figure;
subplot(131); imagesc(x, y, psf);
axis image; title('PSF');
subplot(132); imagesc(x_up, y_up, psf_up);
axis image; title('Upsampled PSF');
subplot(133); imagesc(x_down, y_down, psf_down);
axis image; title('Downsampled PSF');
colormap('hot');

end

function psf = calculate_psf(x, y, p)
% PSF consists of an elliptic 2D Gaussian

% p(1) - amplitude
% p(2) - center x
% p(3) - center y
% p(4) - Standard deviation
% p(5) - constant background
assert(nargin == 3);

sx = p(4) - 0.2;
sy = p(4) + 0.2;

arg_ex = exp(-1/2*((x-p(2))/sx).^2-1/2*((y-p(3))/sy).^2);

psf = p(1) .* arg_ex + p(5); % scale with amplitude and background

end