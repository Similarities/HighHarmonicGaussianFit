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


class FWHM_Gaussian_Fit_Processing:

    def __init__(self, filename, lambda_fundamental, pixel_range, harmonic_selected, filedescription):
        self.filename = filename
        self.ymin = 0
        self.ymax = 2048
        self.xmin = 0
        self.xmax = 2048
        self.pixel_range = pixel_range
        self.picture = np.empty([])
        self.x_backsubstracted = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        self.full_divergence = 17.5
        spatial_detector_calibration = self.full_divergence / 2048
        self.lineout_x = np.arange(self.xmin, self.xmax) * spatial_detector_calibration - self.full_divergence / 2
        self.lineout_y = np.zeros([2048, 1])
        self.filedescription = filedescription
        self.harmonic_selected = harmonic_selected
        self.sigma_temp = float
        self.amplitude_temp = float
        self.center_temp = float
        self.px_boarder = 0
        self.selected_harmonic_in_px()

    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1780:1948], axis=1)
        for x in range(0, self.ymax):
            self.x_backsubstracted[::, x] = self.picture[::, x] - back_mean[x]
        plt.figure(1)
        plt.imshow(self.x_backsubstracted)
        return self.x_backsubstracted

    def grating_function(self):
        x_axis_in_nm = np.empty([2048, 1])
        for x in range(0, self.ymax):
            x_axis_in_nm[x] = 1.24679344e-06 * x ** 2 - 1.65566701e-02 * x + 5.22598053e+01
        return x_axis_in_nm

    def selected_harmonic_in_px(self):
        harmonic_in_nm = self.lambda_fundamental / self.harmonic_selected
        # this function should be inverse of the grating function
        self.px_boarder = (4.71439193e-01 * harmonic_in_nm ** 2 - 1.06651902e+02 * harmonic_in_nm + 4.29603367e+03)
        print(self.px_boarder, "selected N in px")
        return self.px_boarder

    def create_sub_array_px_range(self):
        border_up = int(self.px_boarder - self.pixel_range / 2)
        border_down = int(self.px_boarder + self.pixel_range / 2)
        plt.figure(1)
        plt.hlines(border_up, xmin=0, xmax=2048, color="m", linewidth=0.5)
        plt.hlines(border_down, xmin=0, xmax=2048, color="w", linewidth=1.)
        return self.x_backsubstracted[border_up: border_down, ::]

    def check_fundamental(self):
        sub_array = self.create_sub_array_px_range()
        line_out_y = sub_array[::,1200]
        print(line_out_y)
        line_out_y_1 = list(range(0,self.pixel_range))
        print(line_out_y_1)
        self.plot_x_y(line_out_y_1, line_out_y, 'lineout_over_harmonic_y', 'px', 'counts', 4)
        maximum_in_y = np.where(np.amax(sub_array[::, 1200]))
        print('maximum px position: {0} in px-range {1}'.format(maximum_in_y, self.pixel_range))

    def sum_over_pixel_range(self):
        self.lineout_y = np.sum(self.create_sub_array_px_range(), axis=0)
        return self.lineout_y

    def set_to_zero_offest(self):
        self.lineout_y[::] = self.lineout_y[::] - np.amin(self.lineout_y[150:1900])
        return self.lineout_y

    def fit_gaussian(self):
        self.set_to_zero_offest()
        mod = GaussianModel()
        pars = mod.guess(self.lineout_y, x=self.lineout_x)
        out = mod.fit(self.lineout_y, pars, x=self.lineout_x)
        print(out.fit_report(min_correl=0.15))
        self.sigma_temp = out.params['sigma'].value
        self.amplitude_temp = out.params['amplitude'].value
        self.center_temp = out.params['center'].value
        print('sigma: {0} for N:{1} = {2:8.2f}nm'
              .format(self.sigma_temp, self.harmonic_selected, self.lambda_fundamental / self.harmonic_selected))

        return self.sigma_temp, self.amplitude_temp, self.center_temp

    def plot_x_y(self, x, y, name, x_label, y_label, plot_number):
        plt.figure(plot_number)
        plt.plot(x, y, label=name)
        plt.xlabel('x_label')
        plt.ylabel('y_label')
        plt.legend()

    def plot_fit_function(self):
        xx = np.linspace(-self.full_divergence / 2, self.full_divergence, 1000)
        yy = np.zeros([len(xx), 1])
        for x in range(0, len(xx)):
            yy[x] = (self.amplitude_temp / (self.sigma_temp * ((2 * math.pi) ** 0.5))) * math.exp(
                (-(xx[x] - self.center_temp) ** 2) / (2 * self.sigma_temp ** 2))

        self.plot_x_y(self.lineout_x, self.lineout_y, self.filedescription, 'mrad', 'counts', 2)
        self.plot_x_y(xx, yy, 'fit_function', 'mrad', 'counts', 2)

    def save_linout(self):
        hey = np.stack((self.lineout_x * self.spatial_detector_calibration, self.lineout_y), axis=1)
        np.savetxt(self.filedescription + "lineout" + ".txt", hey, delimiter=' ', fmt='%1.4e')


# insert the following ('filepath/picture_name.tif', fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")

Picture1 = FWHM_Gaussian_Fit_Processing('rotated/spectro1__Wed Jan 23 2019_14.51.31_17.tif', 805., 60, 23, "20190123_17")
Picture1.open_file()
Picture1.background()
Picture1.check_fundamental()

Picture1.sum_over_pixel_range()

Picture1.fit_gaussian()
Picture1.plot_fit_function()
# Picture1.save_data()

plt.show()
