function signal = signal_d( D, ...
                            mf, ...
                            Ib, ...
                            Iu, ...
                            x_bleach, ...
                            y_bleach, ...
                            r_bleach, ...
                            delta_t, ...
                            number_of_pixels, ...
                            number_of_images, ...
                            number_of_pad_pixels)

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
C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 <= (upsampling_factor * r_bleach)^2 ) = Ib;
C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 > (upsampling_factor * r_bleach)^2 ) = Iu;

C0 = imresize(C0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

% FFT of initial condition.
F_C0 = fft2(C0);

% Storage of FFT solution and final solution
F_C = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

signal = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

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
    signal(:, :, t) = abs(ifft2(F_C(:, :, t)));
end
signal = signal(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);

% Take (im)mobile fraction into account and add the free and bound
% contribution to the fluorescence. Even though U0 and B0 are 2-dim
% matrices, they are actually added to each slice in dim 3 of the 3-D
% matrix by writing like this.
signal = mf * signal + (1 - mf) * C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels);

end

