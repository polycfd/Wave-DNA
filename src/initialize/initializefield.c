#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to initialize scalar fields and old scalar fields, where the
latter include the fields for the preceding time steps, needed to the
temporal discretization.
---------------------------------------------------------**/

int InitializeConstScalarField(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *ScalarField, DNA_FLOAT val)
{
  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    ScalarField[iPoint] = val;
  }

  return 0;
}

int InitializeConstOldScalarFields(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT **OldScalarFields, DNA_FLOAT val)
{
  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    OldScalarFields[0][iPoint] = val;
    OldScalarFields[1][iPoint] = val;
  }

  return 0;
}
