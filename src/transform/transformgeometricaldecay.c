#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Function to compute the geometrical factor in the gradient contribution to the Laplacian
for spherical waves. The function is called in <transformeqn.c>
---------------------------------------------------------**/

DNA_FLOAT TransformSphericalDecay(DNA_FLOAT x) { return (2.0 / x); }

DNA_FLOAT TransformNoGeometricalDecay(DNA_FLOAT x) { return (0.0); }