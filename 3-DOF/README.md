Next is a description of the files included in this folder, where the file theodorsen_computations has been implemented or an airfoil and the rest of the files for a finite wing

# pk_strip_theory_aileron_full_span_numerical_mode_shapes 
This script implements a pk method on a 3dof finite wing with a trailing edge aileron along the entire span. It is a limiting case that cannot be performed with the other script included here. The displacements for the mode shapes are numerically integrted t allow for easy manipulations of the geometry.
 

# pk_strip_theory_aileron_not_full_span_numerical_mode_shapes 
This file is very similar to the one described above but allows for the modification of the size of the trailing edge aileron, allowing its size to vary by specifying its initial and final coordinates at the trailing edge, as a percentage of the semispan.

# theodorsen_computations
In this script, the methodolgy implemented by Theodorsen in his paper from 1940: "Mechanism of Flutter A Theoretical and Experimental Investigation of the Flutter Problem" (https://ntrs.nasa.gov/api/citations/19930091762/downloads/19930091762.pdf) to compute flutter has been included.
It enables computation of the flutter point of an airfoil described by non-dimensional parameters. It is an exact procedure to compute the flutter point, providing plots similar to the ones presented in 1940.
