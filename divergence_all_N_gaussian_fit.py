#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019

@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import math

from lmfit.models import GaussianModel


class GaussianFitHighHarmonicDivergence:

    def __init__(self, filename, lambda_fundamental, pixel_range_y, harmonic_selected, file_description):
        self.filename = filename
        # px size full picture * usually 0 - 2048
        self.ymin = 0
        self.ymax = 2048
        # defined Roi in x to avoid boarder effects
        self.xmin = 0
        self.xmax = 2048
        # integration ROI y for each HHG line
        self.ROI_y = pixel_range_y

        self.picture = np.empty([])

        self.x_backsubstracted = np.empty([2048, 2048])

        self.lambdaL = lambda_fundamental

        # calibration of picture in x [full angle], is given with offset here (0 in the middle)
        self.full_divergence = 17.5
        self.C = 17.5 / 2048
        self.lineout_x = np.arange(self.xmin, self.xmax) * self.C - self.full_divergence / 2

        self.lineout_y = np.zeros([2048, 1])

        self.filedescription = file_description

        # defines first harmonic N in pixels, note: the quadratic calibration is not valid for N<10
        self.px_boarder = 0

        self.N_select = harmonic_selected
        self.selected_N_in_px(self.N_select)
        self.N_temp = 0

        self.sigma_temp = float
        self.amplitude_temp = float
        self.center_temp = float

        self.gaussian_result = np.zeros([40, 2])

    def open_file(self):

        self.picture = plt.imread(self.filename)

        return self.picture

    def background(self):

        back_mean = np.mean(self.picture[:, 1780:1948], axis=1)

        i = 1

        N = len(self.picture) - 1

        while i <= N:
            self.x_backsubstracted[::, i] = self.picture[::, i] - back_mean[i]

            i = i + 1

        plt.figure(1)
        plt.xlabel('[px]')
        plt.ylabel('[px]')
        plt.ylim(1000, 650)
        plt.imshow(self.x_backsubstracted)
        plt.legend()

    def selected_N_in_px(self, N):
        self.N_temp = N
        aa = self.lambdaL / self.N_temp

        print(self.lambdaL, "fundamental wavelength, selected harmonic number N:", N, )
        print('= ', aa)

        # this function should be inverse of the grating function
        self.px_boarder = (4.71439193e-01 * (aa ** 2) - 1.06651902e+02 * aa + 4.29603367e+03)

        print(self.px_boarder, "selected N in px")

        # optional tests
        a2 = 1.24679344e-06 * self.px_boarder ** 2 - 1.65566701e-02 * self.px_boarder + 5.22598053e+01
        print("selected N in nm backwards", a2, 'which is N:', round(self.lambdaL / a2))

        return self.px_boarder

    def N_selected_ROI(self):
        # self.selected_N_in_px()

        border_up = int(self.px_boarder - self.ROI_y / 2)
        border_down = int(self.px_boarder + self.ROI_y / 2)

        print(border_up, border_down, "ROI in px")

        subarray = self.x_backsubstracted[border_up: border_down, ::]
        self.lineout_y = np.sum(subarray, axis=0)

        plt.figure(1)

        plt.hlines(border_up, xmin=0, xmax=2048, color="m", linewidth=0.5)

        plt.hlines(border_down, xmin=0, xmax=2048, color="w", linewidth=1.)

        return self.lineout_y

    def fit_gaussian(self):

        mod = GaussianModel()

        pars = mod.guess(self.lineout_y, x=self.lineout_x)
        out = mod.fit(self.lineout_y, pars, x=self.lineout_x)

        print(out.fit_report(min_correl=0.15))

        self.sigma_temp = out.params['sigma'].value
        self.amplitude_temp = out.params['amplitude'].value
        self.center_temp = out.params['center'].value

        plt.figure(2)
        plt.plot(self.lineout_x, self.lineout_y, label=self.filedescription)

        return self.sigma_temp, self.amplitude_temp, self.center_temp

    def plot_fit_function(self):
        xx = np.linspace(-self.full_divergence / 2, self.full_divergence, 1000)
        yy = np.zeros([len(xx), 1])

        for x in range(0, len(xx)):
            yy[x] = (self.amplitude_temp / (self.sigma_temp * ((2 * math.pi) ** 0.5))) * math.exp(
                (-(xx[x] - self.center_temp) ** 2) / (2 * self.sigma_temp ** 2))

        name = self.N_temp

        plt.figure(2)
        plt.plot(xx, yy, label=name, color='c')
        plt.xlabel('[mrad]')
        plt.ylabel('integrated counts')
        plt.legend()

    def batch_over_N(self):

        for x in range(self.N_select, 32):
            self.gaussian_result[x, 0] = x

            self.selected_N_in_px(x)
            self.N_selected_ROI()
            self.fit_gaussian()

            self.plot_fit_function()

            self.gaussian_result[x, 1] = self.sigma_temp

        # clean for empty entries
        self.gaussian_result = np.delete(self.gaussian_result, np.where(~self.gaussian_result.any(axis=1))[0], axis=0)

        plt.figure(3)
        plt.scatter(self.gaussian_result[::, 0], self.gaussian_result[::, 1], label=self.filedescription)
        plt.xlabel('N')
        plt.ylabel('sigma in mrad')
        plt.legend()

        self.save_data()

        return self.gaussian_result

    def save_data(self):

        np.savetxt(self.filedescription + "_div_gaussian" + ".txt", self.gaussian_result, delimiter=' ', fmt='%1.4e')


# insert the following ('filepath/picture_name.tif', fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")

Picture1 = GaussianFitHighHarmonicDivergence('rotated/spectro1__Wed Jan 23 2019_14.49.17_16.tif', 805., 50, 16, "20190123_16")
Picture1.open_file()
Picture1.background()

Picture1.N_selected_ROI()

Picture1.fit_gaussian()
Picture1.plot_fit_function()

Picture1.batch_over_N()

# Picture1.save_data()
