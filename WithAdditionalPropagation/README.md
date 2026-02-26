These MATLAB files implement the original version of the PESCADO method applied to a hydrogen atom exposed to a linearly laser pulse of finite duration within the dipole approximation.
See [https://journals.aps.org/pra/abstract/10.1103/PhysRevA.111.033116](url) for details.

In case the option of calculating the doubly differential ionization probability is selected, differential in both energy and ejection angle, that is, this method involves a 2D interpolation scheeme.
This is done in order to interpolate between the energies pertaining to different L-values. (These pseudo energies are quantized differently for different L-values.) 
