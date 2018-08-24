% function data = signal_db(sys_param, exp_sim_param)

%% Initialization.

% clear
clc
close all hidden

%% Experimental and simulation parameters.

exp_sim_param = struct();

exp_sim_param.pixel_size = 7.5e-07; % m
exp_sim_param.number_of_pixels = 256;

exp_sim_param.number_of_prebleach_frames = 1;
exp_sim_param.number_of_bleach_frames = 2;
exp_sim_param.number_of_postbleach_frames = 3;
exp_sim_param.delta_t = 0.2; % s

exp_sim_param.number_of_pad_pixels = 128;

exp_sim_param.bleach_region.shape = "circular";
exp_sim_param.bleach_region.x = 128; % pixels
exp_sim_param.bleach_region.y = 128; % pixels 
exp_sim_param.bleach_region.r = 15e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.lx = 32%20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.ly = 32%20e-6 / exp_sim_param.pixel_size; % pixels
exp_sim_param.bleach_region.upsampling_factor  = 3;

%% System parameters.

D_SI = 5e-11; % m^2/s
D = D_SI / exp_sim_param.pixel_size^2; % pixels^2 / s
k_on = 1;
k_off = 1;
mobile_fraction = 0.5; % dimensionless
C0 = 1.0; % a.u. original concentration
alpha = 0.6; % a.u.  bleach factor
beta = 1.0; % a.u. imaging bleach factor

sys_param = [D, mobile_fraction, C0, alpha, beta];

% Marginal probabilities of the states.
pi_u = k_off / ( k_on + k_off );
pi_b = k_on / ( k_on + k_off );

%% Simulaton.

tic

% Create a high-resolution bleach mask which is then downsampled to more
% accurately represent edges of the bleach region.
bleach_mask = create_bleach_mask(alpha, exp_sim_param);

% Create imaging bleach mask.
imaging_bleach_mask = create_imaging_bleach_mask(beta, exp_sim_param);

% Fourier space grid squared magnitude.
XSISQ = create_fourier_grid(exp_sim_param);

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

% Prebleach.
U0 = pi_u * C0;
B0 = pi_b * C0;

U_mobile = mobile_fraction * U0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
B_mobile = mobile_fraction * B0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
C_immobile = (1 - mobile_fraction) * C0 * ones(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);

C_prebleach = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_prebleach_frames);

for current_frame = 1:exp_sim_param.number_of_prebleach_frames
    F_U_mobile = fft2(U_mobile);
    F_B_mobile = fft2(B_mobile);
    
    CONST11 = QQ11 .* (QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile);
    CONST12 = QQ12 .* (QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile);
    CONST21 = QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile;
    CONST22 = QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile;
    CONST1 = exp(LL11 * exp_sim_param.delta_t);
    CONST2 = exp(LL22 * exp_sim_param.delta_t);
    
    F_U_mobile = CONST11 .* CONST1 + CONST12 .* CONST2;
    F_B_mobile = CONST21 .* CONST1 + CONST22 .* CONST2;
    
    U_mobile = abs(ifft2(F_U_mobile));
    B_mobile = abs(ifft2(F_B_mobile));
    
    U_mobile = U_mobile .* imaging_bleach_mask;
    B_mobile = B_mobile .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* imaging_bleach_mask;
    
    C_prebleach(:, :, current_frame) = U_mobile + B_mobile + C_immobile;
end

C_prebleach = C_prebleach(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :);

% Bleach.
for current_frame = 1:exp_sim_param.number_of_bleach_frames
    F_U_mobile = fft2(U_mobile);
    F_B_mobile = fft2(B_mobile);
    
    CONST11 = QQ11 .* (QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile);
    CONST12 = QQ12 .* (QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile);
    CONST21 = QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile;
    CONST22 = QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile;
    CONST1 = exp(LL11 * exp_sim_param.delta_t);
    CONST2 = exp(LL22 * exp_sim_param.delta_t);
    
    F_U_mobile = CONST11 .* CONST1 + CONST12 .* CONST2;
    F_B_mobile = CONST21 .* CONST1 + CONST22 .* CONST2;
        
    U_mobile = abs(ifft2(F_U_mobile));
    B_mobile = abs(ifft2(F_B_mobile));
    
    U_mobile = U_mobile .* bleach_mask .* imaging_bleach_mask;
    B_mobile = B_mobile .* bleach_mask .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* bleach_mask .* imaging_bleach_mask;
    
end

% Postbleach.
C_postbleach = zeros(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_postbleach_frames);

for current_frame = 1:exp_sim_param.number_of_postbleach_frames
    F_U_mobile = fft2(U_mobile);
    F_B_mobile = fft2(B_mobile);
    
    CONST11 = QQ11 .* (QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile);
    CONST12 = QQ12 .* (QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile);
    CONST21 = QQinv11 .* F_U_mobile + QQinv12 .* F_B_mobile;
    CONST22 = QQinv21 .* F_U_mobile + QQinv22 .* F_B_mobile;
    CONST1 = exp(LL11 * exp_sim_param.delta_t);
    CONST2 = exp(LL22 * exp_sim_param.delta_t);
    
    F_U_mobile = CONST11 .* CONST1 + CONST12 .* CONST2;
    F_B_mobile = CONST21 .* CONST1 + CONST22 .* CONST2;
    
    U_mobile = abs(ifft2(F_U_mobile));
    B_mobile = abs(ifft2(F_B_mobile));
    
    U_mobile = U_mobile .* imaging_bleach_mask;
    B_mobile = B_mobile .* imaging_bleach_mask;
    
    C_immobile = C_immobile .* imaging_bleach_mask;
    
    C_postbleach(:, :, current_frame) = U_mobile + B_mobile + C_immobile;
    
end

C_postbleach = C_postbleach(exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, exp_sim_param.number_of_pad_pixels + 1:end - exp_sim_param.number_of_pad_pixels, :);

 toc           

% for current_frame = 1:exp_sim_param.number_of_prebleach_frames
%     figure, imagesc([C_prebleach(:, :, current_frame), C_prebleach_sim(:, :, current_frame), C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)]), axis 'equal'
% end
% for current_frame = 1:exp_sim_param.number_of_postbleach_frames
%     figure, imagesc([C_postbleach(:, :, current_frame), C_postbleach_sim(:, :, current_frame), C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)]), axis 'equal'
% end
for current_frame = 1:exp_sim_param.number_of_prebleach_frames
    figure, imagesc(C_prebleach(:, :, current_frame) - C_prebleach_sim(:, :, current_frame)), axis 'equal'
end
for current_frame = 1:exp_sim_param.number_of_postbleach_frames
    figure, imagesc(C_postbleach(:, :, current_frame) - C_postbleach_sim(:, :, current_frame)), axis 'equal'
end
max(abs(C_prebleach(:)-C_prebleach_sim(:)))
max(abs(C_postbleach(:)-C_postbleach_sim(:)))

% figure, hold on, plot(C_postbleach(:, 128), 'k.-'), plot(C_postbleach_sim(:, 128), 'r.-')
% figure, plot(C_postbleach(128, :) - C_postbleach_sim(128, :))
