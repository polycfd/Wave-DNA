#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions construct the central difference approximations of the
first (FDGradient) and the second (FDLaplacian) partial time derivative. As
divergence and gradient only distinguished by context in quasi one-dimensional
coordinates, it is noted that in the present modeling framework, that the first
spatial derivatives to be discretized in the present modeling framework indeed
represent gradients.
---------------------------------------------------------**/

int FDGradient(struct DNA_NumericsFD *NumericsFD, struct DNA_ScalarField *Y, struct DNA_ScalarField *gradField)
{
  int iPoint;
  int endID = NumericsFD->NPoints - 1;

  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adx_order2.dx1[0];
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adx_order2.dx1[1];
  DNA_FLOAT aplus1 = NumericsFD->FDCoeffs.adx_order2.dx1[2];

  // Compute gradient in internal field:
  for (iPoint = 1; iPoint < (NumericsFD->NPoints - 1); iPoint++)
  {
    gradField->val[iPoint] = (amin1 * Y->val[iPoint - 1] + a0 * Y->val[iPoint] + aplus1 * Y->val[iPoint + 1]) / NumericsFD->dXI;
  }

  // Compute gradients at left and right domain boundaries
  gradField->val[0] = (Y->val[1] - Y->val[0]) / NumericsFD->dXI;
  gradField->val[endID] = (Y->val[endID] - Y->val[endID - 1]) / NumericsFD->dXI;

  return 0;
}

int FDLaplacian(struct DNA_NumericsFD *NumericsFD, struct DNA_ScalarField *Y, struct DNA_ScalarField *LaplacianField)
{
  int iPoint;
  int endID = NumericsFD->NPoints - 1;
  DNA_FLOAT dxPow2 = DNA_POW2(NumericsFD->dXI);

  DNA_FLOAT amin1 = NumericsFD->FDCoeffs.adx_order2.dx2[0];
  DNA_FLOAT a0 = NumericsFD->FDCoeffs.adx_order2.dx2[1];
  DNA_FLOAT aplus1 = NumericsFD->FDCoeffs.adx_order2.dx2[2];

  // Compute Laplacian in internal field:
  for (iPoint = 1; iPoint < (NumericsFD->NPoints - 1); iPoint++)
  {
    LaplacianField->val[iPoint] = (amin1 * Y->val[iPoint - 1] + a0 * Y->val[iPoint] + aplus1 * Y->val[iPoint + 1]) / dxPow2;
  }

  // The Laplacian is not required in the BCs but arbitrarily set to zero:
  LaplacianField->val[endID] = 0.0;

  return 0;
}
