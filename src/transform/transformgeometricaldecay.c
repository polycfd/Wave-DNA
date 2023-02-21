#include "DNA.h"
#include "DNA-functions.h"



double TransformSphericalDecay(DNA_FLOAT x)
{
  return (DNA_FLOAT)2.0/x;
}

double TransformNoGeometricalDecay(DNA_FLOAT x)
{
  return (DNA_FLOAT)0.0;
}