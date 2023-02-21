#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to create the grids in the physical domain (GridPhysicalDomain) and in
the computational domain (GridComputationalDomain). The spatial step size is
uniform in both domains. The computational domain ranges from 0 to 1 and it is
temporarily constant. Therefore, it is computed only once upon initialization
(see initializesimulation.c). The physical domain has a time-dependent moving
boundary (xmov) and needs to recomputed at every time step. The number of grid
points (NPoints) is the same in both domains and the time-dependent coordinate 
transformation conveys a one-to-one mapping between the grid points of identical
IDs (iPoint).
---------------------------------------------------------**/

int GridPhysicalDomain(struct DNA_NumericsFD *NumericsFD, struct DNA_Grid *Grid)
{
  int iPoint;
  NumericsFD->dx = (Grid->xfix - Grid->xmov)/((DNA_FLOAT)NumericsFD->NPoints - 1.0);
  
  for(iPoint=0; iPoint<(NumericsFD->NPoints); iPoint++)
  {
    Grid->x[iPoint] = Grid->xmov + ((DNA_FLOAT)iPoint)*(NumericsFD->dx);
  }
  
  return 0;
}


int GridComputationalDomain(struct DNA_NumericsFD *NumericsFD, struct DNA_Grid *Grid)
{
  int iPoint;
  NumericsFD->dXI = 1.0/((DNA_FLOAT) NumericsFD->NPoints - 1.0);
  
  for(iPoint=0; iPoint<(NumericsFD->NPoints); iPoint++)
  {
    Grid->XI[iPoint] = ((DNA_FLOAT)iPoint)*(NumericsFD->dXI);
  }
  
  return 0;
}
