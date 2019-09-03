function example_2d_shift()

%% psf size
size_x = 20;
size_y = 30;

%% derived values
x = single(0 : size_x - 1)';
y = single(0 : size_y - 1);

x_shifted = x - 1.3;
y_shifted = y + 2.7;

%% PSF parameters
psf_parameters = single([100, (size_x-1)/2, (size_y-1)/2, 2, 10]);

%% calculate PSF 
psf = calculate_psf(x, y, psf_parameters);

%% calculate spline coefficients
coefficients = spline_coefficients(psf);

%% generate upsampled psf
psf_shifted = spline_values(coefficients, x_shifted, y_shifted);

%% figure
figure;
subplot(121); imagesc(x, y, psf);
axis image; title('PSF');
subplot(122); imagesc(x_shifted, y_shifted, psf_shifted);
axis image; title('shifted PSF');
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