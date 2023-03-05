#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Function to compute the geometrical factor in the gradient contribution to the Laplacian
for spherical waves. The function is called in <transformeqn.c>
---------------------------------------------------------**/

double TransformSphericalDecay(DNA_FLOAT x)
{
  return (DNA_FLOAT)2.0/x;
}

double TransformNoGeometricalDecay(DNA_FLOAT x)
{
  return (DNA_FLOAT)0.0;
}