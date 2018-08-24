function data = signal_d( ...
    D, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels)

% Bleach region parameters. If param_bleach contains 3 values it is a
% circular bleach region, if 4 it is a rectangular.
x_bleach = param_bleach(1);
y_bleach = param_bleach(2);
if numel(param_bleach) == 3 % Circular.
    r_bleach = param_bleach(3);
else % Rectangular.
    lx_bleach = param_bleach(3);
    ly_bleach = param_bleach(4);
end

% Initial condition. Create a high resolution initial condition which is 
% then downsampled to avoid too sharp edges. Distribute bound and free 
% according to their marginal probabilities.

upsampling_factor = 3;

[X, Y] = meshgrid(1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels), 1:upsampling_factor*(number_of_pixels + 2 * number_of_pad_pixels));
X = X - 0.5;
Y = Y - 0.5;
x_bleach = number_of_pad_pixels + x_bleach;
y_bleach = number_of_pad_pixels + y_bleach;

C0 = zeros(size(X));
if numel(param_bleach) == 3 % Circular.
    C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 <= (upsampling_factor * r_bleach)^2 ) = Ib;
    C0( C0 == 0 ) = Iu;
else % Rectangular.
    C0( X >= upsampling_factor * (x_bleach - 0.5 * lx_bleach) & X <= upsampling_factor * (x_bleach + 0.5 * lx_bleach) & Y >= upsampling_factor * (y_bleach - 0.5 * ly_bleach) & Y <= upsampling_factor * (y_bleach + 0.5 * ly_bleach) ) = Ib;
    C0( C0 == 0 ) = Iu;
end
C0 = imresize(C0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% FFT of initial condition.
F_C0 = fft2(C0);

% Storage of FFT solution and final solution
F_C = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

data = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

% Fourier space grid and squared magnitude, correctly shifted.
[XSI1, XSI2] = meshgrid(-(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1, ...
                        -(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1);
XSI1 = XSI1 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSI2 = XSI2 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSISQ = XSI1.^2 + XSI2.^2;
XSISQ = ifftshift(XSISQ);

% Time evolution in Fourier space.
for t = 1:number_of_images
    T = t * delta_t;
    F_C(:, :, t) = exp( - D * XSISQ * T ) .* F_C0;
end

% Inverse transform.
for t = 1:number_of_images
    data(:, :, t) = abs(ifft2(F_C(:, :, t)));
end
data = data(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);

% Take (im)mobile fraction into account and add the free and bound
% contribution to the fluorescence.
data = mf * data + (1 - mf) * repmat(C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels), [1, 1, number_of_images]);

end