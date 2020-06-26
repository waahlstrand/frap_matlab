import PIL
import numpy as np
from scipy.fftpack import fft2, ifft2
from scipy.ndimage import gaussian_filter

def bleach(c_mobile, c_immobile, X, D, dt, masks, n_frames):

    # Ensure two dimensional
    c_mobile, c_immobile = c_mobile.squeeze(), c_immobile.squeeze()

    C_mobile    = np.zeros((*c_mobile.shape, n_frames))
    C_immobile  = np.zeros((*c_immobile.shape, n_frames))

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

    c_hat = fft2(c)
    c_hat = np.exp(-D * X * dt) * c_hat
    c     = np.abs(ifft2(c_hat))

    return c


def create_imaging_bleach_mask(beta, n_pixels, n_pad_pixels):

    mask = np.ones((n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels))
    mask[n_pad_pixels:-n_pad_pixels, n_pad_pixels:-n_pad_pixels] = beta

    return mask


def create_bleach_mask(alpha, gamma, n_pixels, n_pad_pixels, bleach_region_shape, x, y, lx=None, ly=None, r=None):

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

    if bleach_region_shape == "circle":
        idx_bleach = (X - x)**2 + (Y - y)**2 <= r ** 2

    elif bleach_region_shape == "rectangle":
        idx_bleach = ( X >= x -0.5 * lx ) & (X <= x +0.5 * lx) & (Y >= y -0.5 * ly) & (Y <= y + 0.5*ly)

    bleach_mask_small[idx_bleach] = alpha

    if gamma > 0:
        n = 4 * np.ceil(2 * upsampling_factor * gamma) + 1
        bleach_mask_small = gaussian_filter(bleach_mask_small, sigma= upsampling_factor * gamma, mode='constant')

    
    bleach_mask_small = resize(bleach_mask_small, size= (ub_x - lb_x + 1, ub_y - lb_y + 1))

    bleach_mask = np.ones((n_pixels + 2 * n_pad_pixels, n_pixels + 2 * n_pad_pixels))
    bleach_mask[n_pad_pixels+lb_x-1:n_pad_pixels+ub_x, n_pad_pixels+lb_y-1:n_pad_pixels+ub_y] = bleach_mask_small

    return bleach_mask



def apply_imaging_bleach(c, masks):

    # Apply all masks sequentially
    for mask in masks:
        c *= mask

    return c

def apply_noise(c, a, b):

    sigma = a + b * c
    c = c + np.sqrt(sigma) * np.random.normal(size=c.shape)

    return c


def resize(x, size: int):
    x = convert(x, 0, 255, np.uint8)
    x = np.array(PIL.Image.fromarray(x).resize(size=size, resample=PIL.Image.BOX))
    x = convert(x, 0, 1, np.float64)
    
    return x



def convert(img, target_type_min, target_type_max, target_type):
    imin = img.min()
    imax = img.max()

    a = (target_type_max - target_type_min) / (imax - imin)
    b = target_type_max - a * imax
    new_img = (a * img + b).astype(target_type)
    return new_img