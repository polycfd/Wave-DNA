#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Function to compute the metrics of the coordinate transformation.

detJacobi: Jacobian of the coordinate transformation (dxi/dx).
divq: Divergence due to the fact that the physical can be stretched/compressed,
      so that the grid points are moving away from/towards each other (d^2xi/dtdx).
      This quantity is time-dependent but spatially constant due to the linear
      coordinate transformation.
dxidetJacobi: Included for completeness but always zero due to the linear
              coordinate transformation. Takes non-zero values if
              x(XI) is nonlinear.
---------------------------------------------------------**/

int GridMotion(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{
  Fields->Grid.detJacobi = Fields->Grid.XI[RunOptions->NumericsFD.NPoints-1]/(Fields->Grid.xfix - MovingBoundary->R);
  Fields->Grid.divq = MovingBoundary->U/(Fields->Grid.xfix - MovingBoundary->R);
  Fields->Grid.dxidetJacobi = 0.0;
  
  int iPoint;
  DNA_FLOAT DomainLength = Fields->Grid.xfix - MovingBoundary->R;
  
  for(iPoint=0; iPoint<RunOptions->NumericsFD.NPoints; iPoint++)
  {
    Fields->Grid.q[iPoint] = (Fields->Grid.XI[iPoint] - Fields->Grid.XI[RunOptions->NumericsFD.NPoints-1])/DomainLength*MovingBoundary->U;

    //The following two functions for dtq are equivalent
    //Fields->Grid.dtq[iPoint] = (Fields->Grid.x[iPoint] - Fields->Grid.xfix)*Fields->Grid.detJacobi*(MovingBoundary->Udot*DomainLength + 2.0*DNA_POW2(MovingBoundary->U))/DNA_POW2(DomainLength);
    Fields->Grid.dtq[iPoint] = 2.0*Fields->Grid.q[iPoint]*Fields->Grid.divq 
                              -  Fields->Grid.detJacobi*(Fields->Grid.XI[RunOptions->NumericsFD.NPoints-1] - Fields->Grid.XI[iPoint])*MovingBoundary->Udot;
  }
  
  return 0;
}


