# pk-2-and-3-dof-strip-theory-time-and-frerqueny-domain
Codes used in the development of my Bachelor Thesis to implement a pk method on simple wings.


This repository contains a set of files to analyze 2 and 3 degree of freedom simple aeroelastic systems. 

These files implement a pk analysis over a simple finite wing, following the approach described in Wright and Cooper. Additionally, a time domain simulation of the system is included, employing the codes developed by F.M.A. Mitrotta, for a GRFP system identification algorithm. 

The original code by F.M.A. Mitrotta included in the file GRFP_2DOF can be found in:


https://github.com/fmamitrotta/learn-aeroelastic-sid/blob/main/notebooks/04_Identification_of_a_2-DOF_Mass-Spring-Damper_System_with_the_RFP_Method.ipynb


Next is a description of the file included:

----- timedomain_gust_pk_GRFP -----

This file performs a time domain simulation on a finite 3d wing with 2 degrees of freedom (pitch and plunge), and additionally implements a pk method to characterize the system for a given speed.

The input to the system is a gust with a certain angle acting on the wing. It is randomly generated with a specified seed for reproducibility. The main parameter adjusted to modify the results is the amplitude of this angle, apart from the characteristics of the structure and aerodynamic model. 

This file generates two output files: SIDFiles.mat and FreqResults.mat 
Both of them are used for system identification, fed to GRFP_2DOF.

----- FreqResults -----

This file contains the data from the pk method computed by pk_3d_2dof_RTJones_aero. It includes frequency and damping estimations for each velocity, providing reference values for the GRFP method, and providing a reference.


----- SIDFiles -----
Contains the main variables from the time domain simulation of the system to be used in system identification techniques. The variables are the necessary input for the modal parameters to be obtained with GRFP_2DOF.

----- GRFP_2DOF -----
This file contains code developed by F.M.A. Mitrotta on system identification. This  file is nearly identical to Notebook 4 on system identification, with a few additions to employ it on an aeroelastic SIMO system.

For a given range of speeds, it returns the full implemented GRFP method, with illustrative graphs for the degrees of freedom considered.

-----pk_strip_theory_aileron_full_span_numerical_mode_shapes -----
This script implements a pk method on a 3dof finite wing with a trailing edge aileron along the entire span. It is a limiting case that cannot be performed with the other script included here. The displacements for the mode shapes are numerically integrted t allow for easy manipulations of the geometry.
 

----- pk_strip_theory_aileron_not_full_span_numerical_mode_shapes -----
This file is very similar to the one described above but allows for the modification of the size of the trailing edge aileron, allowing its size to vary by specifying its initial and final coordinates at the trailing edge, as a percentage of the semispan.
