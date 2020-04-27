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








class Open_and_Plot_Picture:
    
    
        def __init__(self, filename,lambdaL,ROI_y,N_select,filedescription):
            self.filename = filename
            # px size full picture * usually 0 - 2048
            self.ymin = 0
            self.ymax = 2048
            # defined Roi in x to avoid boarder effects
            self.xmin = 0
            self.xmax = 2048
            # integration ROI y for each HHG line
            self.ROI_y = ROI_y
            
            self.picture = np.empty([])

            self.x_backsubstracted=np.empty([2048, 2048])
            self.x_axis_in_nm =np.empty([2048,1])

            self.lambdaL = lambdaL



            self.FWHM_for_N = np.zeros([40,3])
            
            #calibration of picture in x [full angle], is given with offset here (0 in the middle)
            self.full_divergence = 17.5
            self.C = 17.5/2048
            self.lineout_x = np.arange(self.xmin,self.xmax)*self.C - self.full_divergence/2
            self.lineout_y = np.zeros([2048,1])


            
            self.filedescription = filedescription
            self.N_select = N_select

            self.px_boarder = 0
            self.selected_N_in_px()
            
            self.sigma_temp = float
            self.amplitude_temp = float
            self.center_temp = float



        def open_file(self):
            
            self.picture = plt.imread(self.filename)
            
            return self.picture
            
      
        
        def background(self):
            
            back_mean=np.mean(self.picture[:, 1780:1948], axis = 1)
            
            
            i=1
            
            N=len(self.picture)-1
            
            while i<= N:
                
                self.x_backsubstracted[::,i] = self.picture[::,i]- back_mean[i]

                i = i+1
                
                
            plt.figure(3)
            plt.ylim(0,500)
            
            plt.imshow(self.x_backsubstracted)



        
        def grating_function(self):
            
            #optional for testing purposes

            N = 2048
            i = 0

            while i <= N-1:
                
                self.x_axis_in_nm[i]= 1.24679344e-06 * i ** 2 - 1.65566701e-02 * i +  5.22598053e+01

                i = i+1
                
            return self.x_axis_in_nm
   
    
    
        def selected_N_in_px(self):
            aa = self.lambdaL/ self.N_select
            
            print(self.lambdaL,  "fundamental wavelength, selected harmonic number N:" ,self.N_select,)
            print('= ', aa)
                


            #this function should be inverse of the grating function
            self.px_boarder = (4.71439193e-01 * aa** 2  - 1.06651902e+02 * aa + 4.29603367e+03)
            print(self.px_boarder, "selected N in px")
            
            #optional tests
            #a2 = 1.24679344e-06 * self.px_boarder ** 2 - 1.65566701e-02 * self.px_boarder +  5.22598053e+01
            #print( "selected N in nm backwards",a2, 'which is N:', round(self.lambdaL/a2))
            
            return self.px_boarder
                
                

                

        def N_selected_ROI(self):
        



            border_up = int(self.px_boarder - self.ROI_y/2)
            border_down = int(self.px_boarder + self.ROI_y/2)
            
            print(border_up, border_down, "ROI in px")
            

            
            
            subarray = self.x_backsubstracted[border_up : border_down, ::]
            self.lineout_y = np.sum(subarray, axis=0)
                
            plt.figure(3)
                
                
            plt.hlines(border_up, xmin = 0, xmax = 2048, color ="m", linewidth =0.5)
                
            plt.hlines(border_down, xmin = 0, xmax = 2048, color ="w", linewidth =1.)
            

            return self.lineout_y

    
    
          
      
        def fit_gaussian(self):
            

            mod = GaussianModel()

            pars = mod.guess(self.lineout_y, x= self.lineout_x)
            out = mod.fit(self.lineout_y, pars, x=self.lineout_x)

            print(out.fit_report(min_correl=0.15))

            self.sigma_temp = out.params['sigma'].value
            self.amplitude_temp = out.params['amplitude'].value
            self.center_temp = out.params['center'].value
            
            
            
            print(self.sigma_temp)
            
            plt.figure(1)
            plt.plot(self.lineout_x, self.lineout_y, label = self.filedescription)
            
            
            return self.sigma_temp, self.amplitude_temp, self.center_temp
        
        
        def plot_fit_function(self):
            xx = np.linspace(-self.full_divergence/2, self.full_divergence, 1000)
            yy = np.zeros([len(xx),1])
            
            for x in range(0,len(xx)): 
                yy[x] = (self.amplitude_temp /(self.sigma_temp*((2*math.pi)**0.5)))* math.exp((-(xx[x]-self.center_temp)**2)/(2*self.sigma_temp**2))

            plt.figure(1)
            plt.plot(xx,yy, label = 'fitfunction')
            plt.legend()
            
            
        
        
        
        def save_data(self):

            hey = np.stack((self.lineout_x*self.C, self.lineout_y), axis=1)

            np.savetxt(self.filedescription+"lineout"+".txt", hey, delimiter=' ', fmt='%1.4e')
  


        


            
     
            






# insert the following ('filepath/picture_name.tif', fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")

Picture1=Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_14.51.31_17.tif', 799., 60, 30,"20190123_17")
Picture1.open_file()
Picture1.background()

Picture1.N_selected_ROI()

Picture1.fit_gaussian()
Picture1.plot_fit_function()
#Picture1.save_data()









