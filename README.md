# Gaussian_fit_HHG_divergence
Gaussian fit for HHG divergence

High Harmonic spectra (HHG) recorded on detector images are evaluated for their spatial divergence property. 
Accordingly the x-Axis corresponds to the spatial axis, where a calibration mrad/px could be implemented. As the HHG divergence is not the same
over teh harmonic number, the evaluation must be done for each harmonic line. Each harmonic line corresponds via a spectral calibration a
certain px position in y-axis on the image, and each harmonic line usually as well contains of a spectral width (here called px_range in y). 
The algorithm provides the following: 

- opens .tif 16 bit pictures in a batch file
- areal background removal (mean value of certain area over the spectral axis (y))
- sum of the signal over px_range in x (ROI)
- sets background to 0 for the line_out
- additional x-depending substraction can be taken
- spectral calibration px to nm (wavelenght) - and vice versa (for nonlinear ROI of integration needed)
- sum of spectral 1/(i+/- something) (i harmonic number N, something < 1) in spectral range over certain number of i
- applies Gaussian fit: used lftm.model - delivering sigma which correspond to the beamwaist (w(z)) for gaussian beams 
for details read wiki about gaussian beams.
- includes certain plots for controlling the results
- saves picture and data in a file 
- can deliver dE/E values for each harmonic number evaluated

Important: 
The gaussian fit can only deliver valid results for a certain minimum dynamic range and
for the case, that the signal distribution is gaussian. 


Python 3.6 

