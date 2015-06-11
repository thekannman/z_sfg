# z_sfg
This program calculates sum-frequency generation (SFG) spectra of aqueous solutions from molecular dynamics simulations. This is done according to the work of the Skinner group. See, for instance, J. Chem. Phys., 2008, 129, 214705.

THIS PROGRAM HAS NOT BEEN THOROUGHLY TESTED!!!

Current limitations include:
* Only trajectories from the GROMACS simulation package can be used as input.

The following libraries are required:
* GROMACS XTC library for reading position/velocity files. 
  http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library
* The Armadillo library for matrix calculations. http://arma.sourceforge.net/
* The Boost libraries for input option parsing, 4D tensor calculations, lexical casting, and output formatting. 
  http://www.boost.org/
