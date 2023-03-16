#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions compute the sum of the acoustic potential Phi times the
finite difference coefficients over previous time steps. This sum is needed
for the explicit finite difference discretization of the temporal derivatives.
It is used in transformeqn.c to construct the terms of the predictor-corrector
method (solvepredictorcorrector.c) based on the transformed govering equation.
It is further used to construct the explicit finite difference approximations
of the first and second partial time derivatives in fddt.c.

Notation: amin is the finite difference coefficient associated with the previous
(old) time level (min1), and amin2 corresponds to the time level preceding the
previous time level (min2).

Remark: In the case of the first-order approximation of the first time derivative
(FDSumFDcoeffTimesPhi_dt1), only the previous time level (min1) constitutes to this
sum. In the case of the first-order approximation of the second time derivative,
there are two time levels involved, i.e.:

d1phi/dt1 = (a0*phinew + amin1*phimin1)/dt
d2phi/dt2 = (a0*phinew + amin1*phimin1 + amin2*phimin2)/dt^2
---------------------------------------------------------**/

int FDSumFDcoeffTimesPhi_dt1(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT **OldPhiFields)
{
  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adt_order1.dt1[1];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    SumFDcoeffTimesPhi[iPoint] = amin1 * OldPhiFields[0][iPoint];
  }

  return 0;
}

int FDSumFDcoeffTimesPhi_dt2(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT **OldPhiFields)
{
  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adt_order1.dt2[1];
  DNA_FLOAT amin2 = NumericsFD->FDCoeffs.adt_order1.dt2[2];

  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    SumFDcoeffTimesPhi[iPoint] = amin1 * OldPhiFields[0][iPoint] + amin2 * OldPhiFields[1][iPoint];
  }

  return 0;
}