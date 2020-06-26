import torch
import numpy as np
from scipy.fft import fftshift, ifftshift


def fourier_grid(n_pixels, n_pad_pixels):

    pi = np.pi

    n_elements = int(np.floor((n_pixels + 2 * n_pad_pixels)/2 + (n_pixels + 2 * n_pad_pixels)/2))

    # Create real grid from linspace 
    x = np.linspace(-(n_pixels + 2 * n_pad_pixels)/2, (n_pixels + 2 * n_pad_pixels)/2 - 1, num=n_elements)
    y = x
    XSI1, XSI2 = np.meshgrid( x, y ) 

    # Normalize
    XSI1 = (2 * pi) / (n_pixels + 2 * n_pad_pixels) * XSI1
    XSI2 = (2 * pi) / (n_pixels + 2 * n_pad_pixels) * XSI2

    # Create Fourier grid
    XSISQ = ifftshift(XSI1 ** 2 + XSI2 ** 2)

    return XSISQ