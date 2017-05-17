function data = signal_db( ...
    D, ...
    k_on, ...
    k_off, ...
    mf, ...
    Ib, ...
    Iu, ...
    param_bleach, ...
    delta_t, ...
    number_of_pixels, ...
    number_of_images, ...
    number_of_pad_pixels)

% Marginal probabilities of the states.
p_u = k_off / ( k_on + k_off );
p_b = k_on / ( k_on + k_off );

% Bleach region parameters. If param_bleach contains 3 values it is a
% circular bleach region, if 4 it is a rectangular.
x_bleach = param_bleach(1);
y_bleach = param_bleach(2);
if numel(param_bleach) == 3 % Circular.
    bleach_region = 'circular';
    r_bleach = param_bleach(3);
else % Rectangular.
    bleach_region = 'rectangular';
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
switch bleach_region 
    case 'circular'
        C0( (X - upsampling_factor * x_bleach).^2 + (Y - upsampling_factor * y_bleach).^2 <= (upsampling_factor * r_bleach)^2 ) = Ib;
        C0( C0 == 0 ) = Iu;
    case 'rectangular'
        C0( X >= upsampling_factor * (x_bleach - 0.5 * lx_bleach) & X <= upsampling_factor * (x_bleach + 0.5 * lx_bleach) & Y >= upsampling_factor * (y_bleach - 0.5 * ly_bleach) & Y <= upsampling_factor * (y_bleach + 0.5 * ly_bleach) ) = Ib;
        C0( C0 == 0 ) = Iu;
end
C0 = imresize(C0, [number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels]);

B0 = p_b * C0;
U0 = p_u * C0;

% FFT of initial conditions.
F_U0 = fft2(U0);
F_B0 = fft2(B0);

% Storage of FFT solution and final solution
F_U = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);
F_B = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

data = zeros(number_of_pixels + 2 * number_of_pad_pixels, number_of_pixels + 2 * number_of_pad_pixels, number_of_images);

% Fourier space grid and squared magnitude, correctly shifted.
[XSI1, XSI2] = meshgrid(-(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1, ...
                        -(number_of_pixels + 2 * number_of_pad_pixels)/2:(number_of_pixels + 2 * number_of_pad_pixels)/2-1);
XSI1 = XSI1 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSI2 = XSI2 * 2 * pi / (number_of_pixels + 2 * number_of_pad_pixels);
XSISQ = XSI1.^2 + XSI2.^2;
XSISQ = ifftshift(XSISQ);

% Precompute elements of the different matrices of the diagonalized form
% of the PDE system matrix, excluding time t which LL will be multiplied 
% with later.
QQ11 = -(k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);
QQ12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*k_on);

LL11 = -((k_on + k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;
LL22 = -((k_on + k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ))./2;

QQinv11 = -k_on./(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2);
QQinv12 = -(k_on - k_off - (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));
QQinv21 = -QQinv11;
QQinv22 = (k_on - k_off + (D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ + k_on^2 + 2.*k_on.*k_off + k_off^2).^(1./2) + D.*XSISQ)./(2.*(2.*k_on.*k_off + k_on^2 + k_off^2 + D^2.*XSISQ.^2 + 2.*D.*k_on.*XSISQ - 2.*D.*k_off.*XSISQ).^(1./2));

% Time evolution in Fourier space.
CONST11 = QQ11 .* (QQinv11 .* F_U0 + QQinv12 .* F_B0);
CONST12 = QQ12 .* (QQinv21 .* F_U0 + QQinv22 .* F_B0);
CONST21 = QQinv11 .* F_U0 + QQinv12 .* F_B0;
CONST22 = QQinv21 .* F_U0 + QQinv22 .* F_B0;

for t = 1:number_of_images
    T = t * delta_t;
    CONST1 = exp(LL11 * T);
    CONST2 = exp(LL22 * T);
    F_U(:, :, t) = CONST11 .* CONST1 + CONST12 .* CONST2;
    F_B(:, :, t) = CONST21 .* CONST1 + CONST22 .* CONST2;
end

% Inverse transform.
for t = 1:number_of_images
    data(:, :, t) = abs(ifft2(F_U(:, :, t) + F_B(:, :, t)));
end
data = data(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels, :);

% Take (im)mobile fraction into account and add the free and bound
% contribution to the fluorescence.
data = mf * data + (1 - mf) * repmat(C0(number_of_pad_pixels+1:end-number_of_pad_pixels, number_of_pad_pixels+1:end-number_of_pad_pixels), [1, 1, number_of_images]);

end