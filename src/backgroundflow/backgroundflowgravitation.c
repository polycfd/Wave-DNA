#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
This function is only applicable in the stationary acoustic black hole sceneraio.
It includes the effect of a gravitational potential on the background flow field.
If not applicable or if not acitvated in the options file, the pointer refers to the dummy.
---------------------------------------------------------**/

double BackgroundFlowGravitationalPotential(struct DNA_RunOptions *RunOptions, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT x)
{
  return (DNA_FLOAT) 2.0 * DNA_POW2(FluidProperties->c0) * RunOptions->radius_horizon / DNA_POW2(x);
}

double BackgroundFlowGravitationalPotential_dummy(struct DNA_RunOptions *RunOptions, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT x)
{
  return (DNA_FLOAT) 0.0;
}