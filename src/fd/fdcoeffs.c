#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Central difference (dx) and forward difference (dt) coefficients. 
---------------------------------------------------------**/

int FDFiniteDifferenceCoeffs(struct DNA_NumericsFD *NumericsFD)
{   
  NumericsFD->FDCoeffs.adt_order1.dt1[0] = 1.0;
  NumericsFD->FDCoeffs.adt_order1.dt1[1] = -1.0;

  NumericsFD->FDCoeffs.adt_order1.dt2[0] = 1.0;
  NumericsFD->FDCoeffs.adt_order1.dt2[1] = -2.0;
  NumericsFD->FDCoeffs.adt_order1.dt2[2] = 1.0;

  NumericsFD->FDCoeffs.adx_order2.dx1[0] = -0.5;
  NumericsFD->FDCoeffs.adx_order2.dx1[1] = 0.0;
  NumericsFD->FDCoeffs.adx_order2.dx1[2] = 0.5;

  NumericsFD->FDCoeffs.adx_order2.dx2[0] = 1.0;
  NumericsFD->FDCoeffs.adx_order2.dx2[1] = -2.0;
  NumericsFD->FDCoeffs.adx_order2.dx2[2] = 1.0;

  return 0;
}