import torchvision
import torch

def bleach(c_mobile, c_immobile, X, D, dt, masks, n_frames):

    # Ensure two dimensional
    c_mobile, c_immobile = c_mobile.squeeze(), c_immobile.squeeze()

    C_mobile    = torch.zeros((*c_mobile.shape, n_frames))
    C_immobile  = torch.zeros((*c_mobile.shape, n_frames))

    # Solve diffusion equation for given number of timesteps
    for frame in range(n_frames):

        # Solve one step in Fourier domain
        c_mobile = fourier_step(c_mobile, X, D, dt)

        # Apply imaging bleach
        c_mobile    = apply_imaging_bleach(c_mobile, masks)
        c_immobile  = apply_imaging_bleach(c_immobile, masks)

        # Save concentration intensities
        C_mobile[:, :, frame]   = c_mobile
        C_immobile[:, :, frame] = c_immobile

    return C_mobile, C_immobile

def fourier_step(c, X, D, dt):

    c_hat = torch.fft(c, signal_ndim=2)
    c_hat = torch.exp(-D * X * dt) * c_hat
    c     = torch.abs(torch.ifft(c_hat, signal_ndim=2))

    return c

def fourier_grid(n_pixels, n_pad_pixels):

    return None
    
def apply_imaging_bleach(c, masks):

    for mask in masks:
        c *= mask

    return c

def apply_noise(c, a, b):

    sigma = a + b * c
    c = c + torch.sqrt(sigma) * torch.randn(c.shape)

    return c
