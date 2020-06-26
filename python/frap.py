import os

import numpy as np

from bleaching import bleach, apply_noise, create_imaging_bleach_mask, create_bleach_mask
from fourier import fourier_grid

class FRAP:

    def __init__(self, **params):
        super(FRAP, self).__init__()

        # System parameters
        #self.D                  = params["system"]["D"]
        #self.c0                 = params["system"]["c0"]
        #self.alpha              = params["system"]["alpha"]
        #self.beta               = params["system"]["beta"]
        #self.gamma              = params["system"]["gamma"]
        #self.mobile_fraction    = params["system"]["mobile_fraction"]
        #self.immobile_fraction  = 1-self.mobile_fraction


        # Experiment parameters
        self.dt                     = params["timestep"]
        self.n_pixels               = params["number_of_pixels"]
        self.n_pad_pixels           = params["number_of_pad_pixels"]
        self.n_prebleach_frames     = params["number_of_prebleach_frames"]
        self.n_bleach_frames        = params["number_of_bleach_frames"]
        self.n_postbleach_frames    = params["number_of_postbleach_frames"]
        self.bleach_region_shape    = params["bleach_region"]["shape"]
        self.center_x               = params["bleach_region"]["x"]
        self.center_y               = params["bleach_region"]["y"]
        self.length_x               = params["bleach_region"]["lx"]
        self.length_y               = params["bleach_region"]["ly"]
        self.radius                 = params["bleach_region"]["r"]
        self.pixel_size             = params["pixel_size"]

        # Scale to pixels
        self.length_x       /= self.pixel_size
        self.length_y       /= self.pixel_size
        self.radius         /= self.pixel_size


        self.X                  = fourier_grid(self.n_pixels, self.n_pad_pixels)
    
    def _create_bleach_mask(self, alpha, gamma):

        return create_bleach_mask(alpha, 
                                  gamma, 
                                  self.n_pixels, self.n_pad_pixels, 
                                  self.bleach_region_shape, 
                                  self.center_x, self.center_y, self.length_x, self.length_y, self.radius)
 
    def _create_imaging_bleach_mask(self, beta):

        return create_imaging_bleach_mask(beta, self.n_pixels, self.n_pad_pixels)

    def _create_bleach_region_indicator(self):

        indicator =  1 - create_bleach_mask(0, 0, 
                                      self.n_pixels, self.n_pad_pixels, 
                                      self.bleach_region_shape, 
                                      self.center_x, self.center_y, self.length_x, self.length_y, self.radius)

        indicator = indicator[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels]
        indicator = indicator / np.sum(indicator.flatten())

        return indicator

    def _apply_noise(self, c, a, b):

        return apply_noise(c, a, b)

    def _bleach(self, C_mobile, C_immobile, D, masks, n_frames):

        return bleach(C_mobile, C_immobile, self.X, D, self.dt, masks, n_frames)


    def _signal(self, **params):
        
        D               = params["D"]
        c0              = params["c0"]
        alpha           = params["alpha"]
        mobile_fraction = params["mobile_fraction"]
        beta            = params["beta"]
        gamma           = params["gamma"]
        a               = params["a"]
        b               = params["b"]

        immobile_fraction = mobile_fraction - 1
        D /= self.pixel_size ** 2

        # Create masks
        imaging_mask = self._create_imaging_bleach_mask(beta)
        bleach_mask  = self._create_bleach_mask(alpha, gamma)

        # Initialize concentrations
        C_mobile_init    = mobile_fraction * c0 * np.ones((self.n_pixels + 2 * self.n_pad_pixels, 
                                                             self.n_pixels + 2 * self.n_pad_pixels))

        C_immobile_init  = immobile_fraction * c0 * np.ones((self.n_pixels + 2 * self.n_pad_pixels, 
                                                               self.n_pixels + 2 * self.n_pad_pixels))


        ######### Prebleach ###########
        C_mobile, C_immobile = self._bleach(C_mobile_init, C_immobile_init, 
                                      D, 
                                      [imaging_mask], 
                                      self.n_prebleach_frames)
        
        # (n_pixels, n_pixels, n_prebleach_frames)
        C_prebleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
                      C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        ########## Bleach #############
        C_mobile, C_immobile = self._bleach(C_mobile[:, :, -1], C_immobile[:, :, -1],
                                      D, 
                                      [bleach_mask, imaging_mask],
                                      self.n_bleach_frames)

        # (n_pixels, n_pixels, n_bleach_frames)
        #C_bleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
        #           C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        ########## Postbleach #############
        C_mobile, C_immobile = self._bleach(C_mobile[:, :, -1], C_immobile[:, :, -1],
                                      D, 
                                      [imaging_mask],
                                      self.n_postbleach_frames)

        # (n_pixels, n_pixels, n_bleach_frames)
        C_postbleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
                       C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        C_prebleach, C_postbleach = self._apply_noise(C_prebleach, a, b), self._apply_noise(C_postbleach, a, b)

        return C_prebleach, C_postbleach

    def recovery_curve(self):

        return None

    def generate(self, mode="pixel", **params):

        C_prebleach, C_postbleach = self._signal(**params)

        if mode == "pixel":
            return C_prebleach, C_postbleach

        elif mode == "rc":
            return self.recovery_curve(C_prebleach, C_postbleach)

        elif mode == "both":
            return C_prebleach, C_postbleach, self.recovery_curve(C_prebleach, C_postbleach)

    def as_generator(self, **params):
        while(True):
            yield self._signal(**params)