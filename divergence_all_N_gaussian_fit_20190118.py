# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019

@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import math

from lmfit.models import GaussianModel


class GaussianFitHighHarmonicDivergence:

    def __init__(self, filename, lambda_fundamental, harmonic_selected):
        self.filename = filename
        # px size full picture * usually 0 - 2048
        self.ymin = 0
        self.ymax = 2048
        # define Roi in x to avoid boarder effects
        self.xmin = 150
        self.xmax = 1450
        # integration ROI y for each HHG line
        self.picture = np.empty([])
        self.picture_background = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        # calibration of picture in x [full angle], is given with offset here (0 in the middle)
        self.full_divergence = 17.5
        self.maximum_harmonic = 29
        self.harmonic_selected = harmonic_selected
        self.lineout_x = self.create_x_axis_in_mrad()
        self.lineout_y = np.zeros([2048, 1])
        self.filedescription = self.filename
        # defines first harmonic N in pixels, note: the quadratic calibration is not valid for N<10
        self.border_up, self.border_down = self.energy_range()
        self.sigma_temp = float
        self.amplitude_temp = float
        self.center_temp = float
        self.gaussian_result = np.zeros([self.maximum_harmonic, 5])

    def create_x_axis_in_mrad(self):
        c = self.full_divergence / 2048
        return np.arange(self.xmin, self.xmax) * c

    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1700:1800], axis=1)
        for x in range(0, 2048):
            self.picture_background[::, x] = self.picture[::, x] - back_mean[x]

        self.background_x()
        plt.figure(1)
        plt.xlabel('[px]')
        plt.ylabel('[px]')
        plt.imshow(self.picture_background, label=self.filedescription)
        return self.picture_background

    def background_x(self):
        back_mean = np.mean(self.picture[1900:2000, :], axis=0)
        for x in range(0, 2048):
            self.picture_background[x, ::] = self.picture[x, ::] - back_mean[x]
        return self.picture_background

    def nm_in_px(self, energy_nm):
        # this function should be inverse of the grating function
        self.px_boarder = np.rint(3.71756520e-01 * energy_nm ** 2 - 9.89219185e+01 * energy_nm + 4.50769013e+03)
        return self.px_boarder

    def energy_range(self):
        print(self.harmonic_selected)
        previous_harmonic = self.lambda_fundamental / (self.harmonic_selected - 0.2)
        next_harmonic = self.lambda_fundamental / (self.harmonic_selected + 0.2)
        print(previous_harmonic, next_harmonic)
        self.border_up = np.int(self.nm_in_px(previous_harmonic))
        self.border_down = np.int(self.nm_in_px(next_harmonic))
        print(self.border_up, self.border_down, "ROI in px")
        self.pixel_range = np.int(self.border_down - self.border_up)
        print(self.pixel_range, 'ROI in pixel range')
        return self.border_up, self.border_down

    def delta_energy(self):
        lower = int(self.border_up)
        upper = int(self.border_down)
        delta = self.px_in_nm(lower) - self.px_in_nm(upper)
        delta_vs_energy = delta / (self.lambda_fundamental / self.harmonic_selected)
        return delta_vs_energy, delta

    def px_in_nm(self, px_number):
        return 1.22447518e-06 * px_number ** 2 - 1.73729829e-02 * px_number + 5.82820234e+01

    def create_sub_array_px_range(self):
        border_up = int(self.border_up)
        border_down = int(self.border_down)
        self.print_h_lines()
        return self.picture_background[border_up: border_down, self.xmin:self.xmax]

    def print_h_lines(self):
        fundamental_in_px = self.nm_in_px(self.lambda_fundamental / self.harmonic_selected)
        plt.figure(1)
        plt.hlines(self.border_up, xmin=0, xmax=2048, color="b", linewidth=0.3)
        plt.hlines(self.border_down, xmin=0, xmax=2048, color="w", linewidth=0.3)
        # plt.hlines(fundamental_in_px, xmin=0, xmax=2048, color='r', linewidth=0.1)
        plt.vlines(self.xmin, ymin=0, ymax=2048, color="w", linewidth=0.5)
        plt.vlines(self.xmax, ymin=0, ymax=2048, color="w", linewidth=0.5)

    def check_fundamental(self):
        sub_array = self.create_sub_array_px_range()
        line_out_y = sub_array[::, 1200]
        line_out_y_1 = np.arange(0, self.ymax)
        self.plot_x_y(line_out_y_1, line_out_y, 'lineout_over_harmonic_y', 'px', 'counts', 4)
        maximum_in_y = np.where(np.amax(sub_array[::, 1200]))
        print('maximum px position: {0} in px-range {1}'.format(maximum_in_y, self.ymax))

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
        self.sigma_temp = out.params['sigma'].value
        self.amplitude_temp = out.params['amplitude'].value
        self.center_temp = out.params['center'].value
        # print('sigma: {0} for N:{1} = {2:8.2f}nm'
        #    .format(self.sigma_temp, self.harmonic_selected, self.lambda_fundamental / self.harmonic_selected))
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
        # todo: is sigma 1/e^2 - which corresponds to w(0) beamwaist = half beam aperture
        c = self.full_divergence / 2048
        xx = np.linspace(self.xmin * c, self.xmax * c, 1000)
        yy = np.zeros([len(xx), 1])
        for x in range(0, len(xx)):
            yy[x] = (self.amplitude_temp / (self.sigma_temp * ((2 * math.pi) ** 0.5))) * math.exp(
                (-(xx[x] - self.center_temp) ** 2) / (2 * self.sigma_temp ** 2))
        self.plot_x_y(self.lineout_x, self.lineout_y, str(self.harmonic_selected), 'mrad', 'counts', 2)
        self.plot_x_y(xx, yy, 'fit_' + str(self.harmonic_selected), 'mrad', 'counts', 2)

    def batch_over_N(self):
        for x in range(self.harmonic_selected, self.maximum_harmonic):
            self.gaussian_result[x, 0] = x
            self.harmonic_selected = x
            self.border_up, self.border_down = self.energy_range()
            self.fit_gaussian()
            # self.plot_fit_function()
            plt.figure(1)
            plt.hlines(self.nm_in_px(self.lambda_fundamental / self.harmonic_selected), xmin=0, xmax=2048,
                       linewidth=0.5, alpha=0.1)
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
        header_names = (['harmonic number', 'mrad', 'integrated counts in delta E', 'harmonic in nm', 'delta E/E'])
        parameter_info = (
            ['fundamental_nm:', str(self.lambda_fundamental), 'pixel_range:', str(self.ymax), 'xxxx'])
        return np.vstack((header_names, self.gaussian_result, parameter_info))

    def save_data(self):
        result = self.prepare_header()
        print('saved data')
        np.savetxt(self.filedescription[31:42] + '_' + self.filedescription[-6:-4] + ".txt", result, delimiter=' ',
                   header='string', comments='',
                   fmt='%s')


def get_file_list(path_picture):
    tif_files = []
    counter = 0
    for file in os.listdir(path_picture):
        print(file)
        try:
            if file.endswith(".tif"):
                tif_files.append(str(file))
                counter = counter + 1
            else:
                print("only other files found")
        except Exception as e:
            raise e
    return tif_files


def process_files(my_files, path):
    for x in range(20, 21):
        file = path + '/' + my_files[x]
        ProcessImage = GaussianFitHighHarmonicDivergence(file, 800,
                                                         22)
        ProcessImage.open_file()
        ProcessImage.background()
        ProcessImage.batch_over_N()
        ProcessImage.save_data()
        plt.show()


# insert the following ('filepath/picture_name.tif', fundamental frequency (float),
# harmonic number (int, first harmonic), "picture name for plot")


path_picture = 'rotated_20190118'
my_files = get_file_list(path_picture)
process_files(my_files, path_picture)
