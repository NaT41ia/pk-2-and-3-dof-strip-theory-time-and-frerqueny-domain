Next is a description of the file included:

# timedomain_gust_pk_GRFP 

This file performs a time domain simulation on a finite 3d wing with 2 degrees of freedom (pitch and plunge), and additionally implements a pk method to characterize the system for a given speed.

The input to the system is a gust with a certain angle acting on the wing. It is randomly generated with a specified seed for reproducibility. The main parameter adjusted to modify the results is the amplitude of this angle, apart from the characteristics of the structure and aerodynamic model. 

This file generates two output files: SIDFiles.mat and FreqResults.mat 
Both of them are used for system identification, fed to GRFP_2DOF.

# FreqResults

This file contains the data from the pk method computed by pk_3d_2dof_RTJones_aero. It includes frequency and damping estimations for each velocity, providing reference values for the GRFP method, and providing a reference.


# SIDFiles
Contains the main variables from the time domain simulation of the system to be used in system identification techniques. The variables are the necessary input for the modal parameters to be obtained with GRFP_2DOF.

# GRFP_2DOF 
This file contains code developed by F.M.A. Mitrotta on system identification. This  file is nearly identical to Notebook 4 on system identification, with a few additions to employ it on an aeroelastic SIMO system.

For a given range of speeds, it returns the full implemented GRFP method, with illustrative graphs for the degrees of freedom considered.

The original code by F.M.A. Mitrotta included in the file GRFP_2DOF can be found in:
https://github.com/fmamitrotta/learn-aeroelastic-sid/blob/main/notebooks/04_Identification_of_a_2-DOF_Mass-Spring-Damper_System_with_the_RFP_Method.ipynb

# pk_3d_2_dof_Thodorsen_aero
This is a pk method implemented for a finite wing. It provides the flutter point of a model similar to that described by Wright and Cooper but implementing Theodorsen's unsteady model.

# pk_3d_2_dof_RTJones_aero
This is a pk method implemented for a finite wing. It provides the flutter point of a model similar to that described by Wright and Cooper but implementing R.T.Jones approximation to Theodorsen's unsteady model.
