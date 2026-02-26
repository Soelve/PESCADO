These MATLAB files implement the "autnomous" version of the PESCADO. 
It is autonomous in the sense that it allows us to calculate photoelectron spectra using a semi-analytical method when the Hamiltonien no longer carries any time-dependence.
See [https://arxiv.org/abs/2510.05776](url) for details.

There are two alternative ways to calculate energy spectra: Either by using the pseudo energies on the grid - with interpolation - or by using the "true" continuum states generated using the Specal Functions in Physics Toolbox for MATLAB:
[https://se.mathworks.com/matlabcentral/fileexchange/88041-special-functions-in-physics-specfunphys-toolbox](url)

The latter has the advantage that it evades the need to perform a two-dimensional interpolation when calculating the doubly differential photoelectron spectra. The energy eigen function must be pre-calculated and stored in a file; 
the name of this file should be given as input.
