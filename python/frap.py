import os

import torch

from .bleaching import bleach, apply_noise, fourier_grid

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
        self.a                      = params["experiment"]["a"]
        self.b                      = params["experiment"]["b"]
        self.dt                     = params["experiment"]["timestep"]
        self.n_pixels               = params["experiment"]["number_of_pixels"]
        self.n_pad_pixels           = params["experiment"]["number_of_pad_pixels"]
        self.n_prebleach_frames     = params["experiment"]["number_of_prebleach_frames"]
        self.n_bleach_frames        = params["experiment"]["number_of_bleach_frames"]
        self.n_postbleach_frames    = params["experiment"]["number_of_postbleach_frames"]

        self.bleach_mask        = self._create_bleach_mask(**params["experiment"])
        self.X                  = fourier_grid(self.n_pixels, self.n_pad_pixels)
    
    def _create_bleach_mask(self):

        return None

    def _create_imaging_mask(self):

        return None

    def _apply_noise(self, c):

        return apply_noise(c, self.a, self.b)

    def _bleach(self, C_mobile, C_immobile, D, masks, n_frames):

        return bleach(C_mobile, C_immobile, self.X, D, self.dt, masks, n_frames)


    def _signal(self, **params):

        # Create masks
        imaging_mask = self._create_imaging_mask(params["alpha"])
        bleach_mask  = self._create_bleach_mask(params["beta"])

        # Initialize concentrations
        C_mobile_init    = params["mobile_fraction"] * params["c0"] * torch.ones(self.n_pixels + 2 * self.n_pad_pixels, 
                                                                       self.n_pixels + 2 * self.n_pad_pixels)

        C_immobile_init  = (params["mobile_fraction"]-1) * params["c0"] * torch.ones(self.n_pixels + 2 * self.n_pad_pixels, 
                                                                       self.n_pixels + 2 * self.n_pad_pixels)

        ######### Prebleach ###########
        C_mobile, C_immobile = self._bleach(C_immobile_init, C_immobile_init, 
                                      params["D"], 
                                      [imaging_mask], 
                                      self.n_prebleach_frames)
        
        # (n_pixels, n_pixels, n_prebleach_frames)
        C_prebleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
                      C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        ########## Bleach #############
        C_mobile, C_immobile = self._bleach(C_mobile[:, :, -1], C_immobile[:, :, -1],
                                      params["D"], 
                                      [bleach_mask, imaging_mask],
                                      self.n_bleach_frames)

        # (n_pixels, n_pixels, n_bleach_frames)
        #C_bleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
        #           C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        ########## Postbleach #############
        C_mobile, C_immobile = self._bleach(C_mobile[:, :, -1], C_immobile[:, :, -1],
                                      params["D"], 
                                      [imaging_mask],
                                      self.n_postbleach_frames)

        # (n_pixels, n_pixels, n_bleach_frames)
        C_postbleach = C_mobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :] + \
                       C_immobile[self.n_pad_pixels:-self.n_pad_pixels, self.n_pad_pixels:-self.n_pad_pixels, :]

        return C_prebleach, C_postbleach

    def recovery_curve(self):

        return None

    def generate(self, mode="pixel", **params):

        C_prebleach, C_postbleach = self._signal(**params)
        C_prebleach, C_postbleach = self._apply_noise(C_prebleach), self._apply_noise(C_postbleach)

        if mode == "pixel":
            return C_prebleach, C_postbleach

        elif mode == "rc":
            return self.recovery_curve(C_prebleach, C_postbleach)

        elif mode == "both":
            return C_prebleach, C_postbleach, self.recovery_curve(C_prebleach, C_postbleach)