#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
This function calls the pointers referring to the boundary conditions.
---------------------------------------------------------**/

int BoundaryConditions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties)
{
  RunOptions->FDBC_East(RunOptions, Fields, FluidProperties, &Fields->phi, &Fields->Old_phi.o[0]);
  RunOptions->FDBC_West(RunOptions, Fields, FluidProperties, &Fields->phi, &Fields->Old_phi.o[0]);

  return 0;
}