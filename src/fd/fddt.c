#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions construct the finite difference approximations of the
first (FDdt1) and the second (FDdt1) partial time derivative based on the sum
over the previous time steps (SumAY). SumAy is computed in fdsumay.c.
---------------------------------------------------------**/

int FDdt1(struct DNA_NumericsFD *NumericsFD, struct DNA_ScalarField *SumAY, struct DNA_ScalarField *Y, struct DNA_ScalarField *dt1Field)
{
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adt_order1.dt1[0];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    dt1Field->val[iPoint] = (a0 * Y->val[iPoint] + SumAY->val[iPoint]) / NumericsFD->dt;
  }

  return 0;
}

int FDdt2(struct DNA_NumericsFD *NumericsFD, struct DNA_ScalarField *SumAY, struct DNA_ScalarField *Y, struct DNA_ScalarField *dt2Field)
{
  DNA_FLOAT dtPow2 = DNA_POW2(NumericsFD->dt);
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adt_order1.dt2[0];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    dt2Field->val[iPoint] = (a0 * Y->val[iPoint] + SumAY->val[iPoint]) / dtPow2;
  }

  return 0;
}