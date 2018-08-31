function [C_prebleach, C_postbleach] = signal_db(sys_param, exp_sim_param)

% Extract system parameters.
D = sys_param(1);
k_on = sys_param(2);
k_off = sys_param(3);
mobile_fraction = sys_param(4);
C0 = sys_param(5);
alpha = sys_param(6);
beta = sys_param(7);
gamma = sys_param(8);

% Marginal probabilities of the states.
pi_u = k_off / ( k_on + k_off );
pi_b = k_on / ( k_on + k_off );

% Create a high-resolution bleach mask which is then downsampled to more
% accurately represent edges of the bleach region.
bleach_mask = create_bleach_mask(alpha, gamma, exp_sim_param);

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

end