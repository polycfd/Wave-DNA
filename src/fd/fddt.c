#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions construct the finite difference approximations of the
first (FDdt1) and the second (FDdt1) partial time derivative based on the sum
over the previous time steps (SumAY). SumAy is computed in fdsumay.c.
---------------------------------------------------------**/

int FDdt1(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumAY, DNA_FLOAT *Y, DNA_FLOAT *dt1Field)
{
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adt_order1.dt1[0];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    dt1Field[iPoint] = (a0 * Y[iPoint] + SumAY[iPoint]) / NumericsFD->dt;
  }

  return 0;
}

int FDdt2(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumAY, DNA_FLOAT *Y, DNA_FLOAT *dt2Field)
{
  DNA_FLOAT dtPow2 = DNA_POW2(NumericsFD->dt);
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adt_order1.dt2[0];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    dt2Field[iPoint] = (a0 * Y[iPoint] + SumAY[iPoint]) / dtPow2;
  }

  return 0;
}