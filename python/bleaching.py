import torchvision
import torch
import numpy as np
import scipy as scp

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

    pi = np.pi

    n_elements = int(np.floor((n_pixels + 2 * n_pad_pixels)/2 + (n_pixels + 2 * n_pad_pixels)/2))

    # Create real grid from linspace 
    x = np.linspace(-(n_pixels + 2 * n_pad_pixels)/2, (n_pixels + 2 * n_pad_pixels)/2 - 1, num=n_elements)
    y = x
    XSI1, XSI2 = np.meshgrid( x, y ) 

    # Normalize
    XSI1 *= (2 * pi) / (n_pixels + 2 * n_pad_pixels)
    XSI2 *= (2 * pi) / (n_pixels + 2 * n_pad_pixels)

    # Create Fourier grid
    XSISQ = np.fft.ifftshift(XSI1 ** 2 + XSI2 ** 2)

    return torch.Tensor(XSISQ)

def create_imaging_bleach_mask(beta, n_pixels, n_pad_pixels):

    mask = torch.ones(n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels)
    mask[n_pad_pixels:-n_pad_pixels, n_pad_pixels:-n_pad_pixels] = beta

    return mask


def create_bleach_mask(alpha, gamma, bleach_region_shape, x, y, lx=None, ly=None, r=None):

    upsampling_factor=15

    if bleach_region_shape == "circle":

        lb_x = int(np.floor(0.5 + x - r - 8 * gamma))
        ub_x = int(np.ceil(0.5 + x + r + 8 * gamma))
        lb_y = int(np.floor(0.5 + y - r - 8 * gamma))
        ub_y = int(np.ceil(0.5 + y + r + 8 * gamma))

    elif bleach_region_shape == "rectangle":

        lb_x = int(np.floor(0.5 + x - 0.5 * lx - 8 * gamma))
        ub_x = int(np.ceil(0.5 + x + 0.5 * lx + 8 * gamma))
        lb_y = int(np.floor(0.5 + y - 0.5 * ly - 8 * gamma))
        ub_y = int(np.ceil(0.5 + y + 0.5 * ly + 8 * gamma))

    else:
        raise ValueError("Bleach region shape must be circle or rectangle.")

    lb_x = lb_x - 1
    ub_x = ub_x + 1
    lb_y = lb_y - 1
    ub_y = ub_y + 1

    xx = np.linspace(lb_x - 1, ub_x, int(upsampling_factor * (ub_x - lb_x + 1)))
    yy = np.linspace(lb_y - 1, ub_y, int(upsampling_factor * (ub_y - lb_y + 1)))

    [X, Y] = np.meshgrid(xx, yy)

    bleach_mask_small = np.ones(X.shape)

    if bleach_region_shape is "circle":
        idx_bleach = (X - x)**2 + (Y - y)**2 <= r ** 2

    elif bleach_region_shape is "rectangle":
        idx_bleach = ( X >= x -0.5 * lx ) & (X <= x +0.5 * lx) & (Y >= y -0.5 * ly) & (Y <= y + 0.5*ly)

    bleach_mask_small[idx_bleach] = alpha

    if gamma > 0:
        n = 4 * np.ceil(2 * upsampling_factor * gamma) + 1
        bleach_mask_small = scp.ndimage.gaussian_filter(bleach_mask_small, sigma=n, mode='constant')




def apply_imaging_bleach(c, masks):

    # Apply all masks sequentially
    for mask in masks:
        c *= mask

    return c

def apply_noise(c, a, b):

    sigma = a + b * c
    c = c + torch.sqrt(sigma) * torch.randn(c.shape)

    return c
