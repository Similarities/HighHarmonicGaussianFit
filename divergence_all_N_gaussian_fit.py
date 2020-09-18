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
        # define Roi in x to avoid boarder effects
        self.xmin = 0
        self.xmax = 2048
        # integration ROI y for each HHG line
        self.pixel_range_y = pixel_range_y
        self.picture = np.empty([])
        self.picture_background = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        # calibration of picture in x [full angle], is given with offset here (0 in the middle)
        self.full_divergence = 17.5
        self.maximum_harmonic = 33
        self.lineout_x = self.create_x_axis_in_mrad()
        self.lineout_y = np.zeros([2048, 1])
        self.filedescription = file_description
        # defines first harmonic N in pixels, note: the quadratic calibration is not valid for N<10
        self.px_boarder = 0
        self.harmonic_selected = harmonic_selected
        self.selected_harmonic_in_px(self.harmonic_selected)
        self.sigma_temp = float
        self.amplitude_temp = float
        self.center_temp = float
        self.gaussian_result = np.zeros([self.maximum_harmonic, 5])

    def create_x_axis_in_mrad(self):
        c = self.full_divergence / self.xmax
        return np.arange(self.xmin, self.xmax) * c - self.full_divergence / 2

    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1780:1848], axis=1)
        for x in range(0, 2048):
            self.picture_background[::, x] = self.picture[::, x] - back_mean[x]
        plt.figure(1)
        plt.xlabel('[px]')
        plt.ylabel('[px]')
        plt.imshow(self.picture_background, label=self.filedescription)
        return self.picture_background

    def selected_harmonic_in_px(self, N):
        harmonic_in_nm = self.lambda_fundamental / self.harmonic_selected
        # this function should be inverse of the grating function
        self.px_boarder = (4.71439193e-01 * harmonic_in_nm ** 2 - 1.06651902e+02 * harmonic_in_nm + 4.29603367e+03)
        return self.px_boarder

    def delta_energy(self):
        lower = int(self.px_boarder - self.pixel_range_y / 2)
        upper = int(self.px_boarder + self.pixel_range_y / 2)
        delta = self.px_in_nm(lower) - self.px_in_nm(upper)
        delta_vs_energy = delta / self.px_in_nm(self.px_boarder)
        return delta_vs_energy, delta

    def px_in_nm(self, px_number):
        return 1.24679344e-06 * px_number ** 2 - 1.65566701e-02 * px_number + 5.22598053e+01

    def create_sub_array_px_range(self):
        border_up = int(self.px_boarder - self.pixel_range_y / 2)
        border_down = int(self.px_boarder + self.pixel_range_y / 2)
        plt.figure(1)
        plt.hlines(border_up, xmin=0, xmax=2048, color="m", linewidth=0.5)
        plt.hlines(border_down, xmin=0, xmax=2048, color="w", linewidth=1.)
        return self.picture[border_up: border_down, ::]

    def check_fundamental(self):
        sub_array = self.create_sub_array_px_range()
        line_out_y = sub_array[::, 1200]
        line_out_y_1 = np.arange(0, self.pixel_range_y)
        self.plot_x_y(line_out_y_1, line_out_y, 'lineout_over_harmonic_y', 'px', 'counts', 4)
        maximum_in_y = np.where(np.amax(sub_array[::, 1200]))
        print('maximum px position: {0} in px-range {1}'.format(maximum_in_y, self.pixel_range_y))

    def sum_over_pixel_range(self):
        self.lineout_y = np.sum(self.create_sub_array_px_range(), axis=0)
        return self.lineout_y

    def set_to_zero_offest(self):
        self.lineout_y[::] = self.lineout_y[::] - np.amin(self.lineout_y)
        return self.lineout_y

    def fit_gaussian(self):
        self.sum_over_pixel_range()
        self.set_to_zero_offest()
        mod = GaussianModel()
        pars = mod.guess(self.lineout_y, x=self.lineout_x)
        out = mod.fit(self.lineout_y, pars, x=self.lineout_x)
        #print(out.fit_report(min_correl=0.15))
        self.sigma_temp = out.params['sigma'].value
        self.amplitude_temp = out.params['amplitude'].value
        self.center_temp = out.params['center'].value
        #print('sigma: {0} for N:{1} = {2:8.2f}nm'
        #    .format(self.sigma_temp, self.harmonic_selected, self.lambda_fundamental / self.harmonic_selected))
        self.plot_x_y(self.lineout_x, self.lineout_y, str(self.harmonic_selected), ' ', ' ', 2)
        self.plot_fit_function()
        return self.sigma_temp, self.amplitude_temp, self.center_temp

    def plot_x_y(self, x, y, name, x_label, y_label, plot_number):
        plt.figure(plot_number)
        plt.plot(x, y, label=name)
        plt.xlabel(str(x_label))
        plt.ylabel(str(y_label))
        plt.legend()

    def plot_scatter(self, x, y, name, x_label, y_label, plot_number):
        plt.figure(plot_number)
        plt.scatter(x, y, label=name)
        plt.xlabel(str(x_label))
        plt.ylabel(str(y_label))
        plt.legend()

    def plot_fit_function(self):
        xx = np.linspace(-self.full_divergence/2, self.full_divergence/2, 1000)
        yy = np.zeros([len(xx), 1])
        for x in range(0, len(xx)):
            yy[x] = (self.amplitude_temp / (self.sigma_temp * ((2 * math.pi) ** 0.5))) * math.exp(
                (-(xx[x] - self.center_temp) ** 2) / (2 * self.sigma_temp ** 2))
        self.plot_x_y(self.lineout_x, self.lineout_y, str(self.harmonic_selected), 'mrad', 'counts', 2)
        self.plot_x_y(xx, yy, 'fit_'+str(self.harmonic_selected), 'mrad', 'counts', 2)

    def batch_over_N(self):
        for x in range(self.harmonic_selected, self.maximum_harmonic):
            self.gaussian_result[x, 0] = x
            self.harmonic_selected = x
            self.selected_harmonic_in_px(self.harmonic_selected)
            self.fit_gaussian()
            #self.plot_fit_function()
            self.gaussian_result[x, 1] = self.sigma_temp
            self.gaussian_result[x, 2] = np.sum(self.lineout_y[::])
            self.gaussian_result[x, 4], self.gaussian_result[x, 3] = self.delta_energy()

        # clean for empty entries
        self.gaussian_result = np.delete(self.gaussian_result, np.where(~self.gaussian_result.any(axis=1))[0], axis=0)
        self.plot_scatter(self.gaussian_result[::, 0], self.gaussian_result[::, 1], self.filedescription,
                          'harmonic number N', 'divergence in mrad', 3)
        return self.gaussian_result

    def prepare_header(self):
        # insert header line and change index
        header_names = (['harmonic number', 'mrad', 'integrated counts in delta E' ,'harmonic in nm' ,'delta E/E'])
        return np.vstack((header_names, self.gaussian_result))


    def save_data(self):
        result = self.prepare_header()
        np.savetxt(self.filedescription + "_div_gaussian" + ".txt", result, delimiter=' ', header='string', comments='', fmt='%s')


# insert the following ('filepath/picture_name.tif', fundamental frequency (float), pixel_range_y (delta energy), harmonic number (int), "picture name for plot")

Picture1 = GaussianFitHighHarmonicDivergence('rotated/spectro1__Wed Jan 23 2019_16.01.29_36.tif', 801., 30, 16,
                                             "20190123_16")
Picture1.open_file()
Picture1.background()
Picture1.batch_over_N()
Picture1.save_data()
plt.show()

# Picture1.save_data()
