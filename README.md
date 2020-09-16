# Gaussian_fit_HHG_divergence
Gaussian fit for HHG divergence

High Harmonic spectras (HHG) recorded on detector images are evaluated for their spatial divergence property. 
Accordingly the x-Axis corresponds to the spatial axis, where a calibration mrad/px could be implemented. As the HHG divergence is not the same
over teh harmonic number, the evaluation must be done for each harmonic line. Each harmonic line corresponds via a spectral calibration a
certain px position in y-axis on the image, and each harmonic line usually as well contains of a spectral width (here called px_range in y). 
The algorithm provides the following: 
opens .tif 16 bit picture //
areal background removal (mean value of certain area over the spectral axis (y)) // test the px_range and fundamental frequency on the image (nm to px calculation) //
summation of the signal over px_range // sets background to 0 for the line_out// applies Gaussian fit: used lftm.model// includes certain plots for controling the results.

Python 3.6 

