#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions compute the sum of some scalar quantity Y times the
finite difference coefficients over previous time steps. This sum is needed
for the explicit finite difference discretization of the temporal derivatives.
It is used in transformeqn.c to construct the terms of the predictor-corrector
method (solvepredictorcorrector.c) based on the transformed govering equation.
It is further used to construct the explicit finite difference approximations
of the first and second partial time derivatives in fddt.c.
---------------------------------------------------------**/

int FDSumAY_dt1(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumAY, DNA_FLOAT **OldScalarFields)
{
  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adt_order1.dt1[1];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    SumAY[iPoint] = amin1 * OldScalarFields[0][iPoint];
  }

  return 0;
}

int FDSumAY_dt2(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumAY, DNA_FLOAT **OldScalarFields)
{
  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adt_order1.dt2[1];
  DNA_FLOAT amin2 = NumericsFD->FDCoeffs.adt_order1.dt2[2];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    SumAY[iPoint] = amin1 * OldScalarFields[0][iPoint] + amin2 * OldScalarFields[1][iPoint];
  }

  return 0;
}