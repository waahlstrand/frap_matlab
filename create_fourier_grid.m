function XSISQ = create_fourier_grid(exp_sim_param)

[XSI1, XSI2] = ndgrid(  -(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels)/2:(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels)/2-1, ...
                        -(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels)/2:(exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels)/2-1);
XSI1 = XSI1 * 2 * pi / (exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
XSI2 = XSI2 * 2 * pi / (exp_sim_param.number_of_pixels + 2 * exp_sim_param.number_of_pad_pixels);
XSISQ = ifftshift(XSI1.^2 + XSI2.^2);

end

