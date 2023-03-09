#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to initialize scalar fields and old scalar fields, where the
latter include the fields for the preceding time steps, needed to the
temporal discretization.
---------------------------------------------------------**/

int InitializeConstScalarField(struct DNA_NumericsFD *NumericsFD, struct DNA_ScalarField *ScalarField, DNA_FLOAT val)
{
  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    ScalarField->val[iPoint] = val;
  }

  return 0;
}

int InitializeConstOldScalarFields(struct DNA_NumericsFD *NumericsFD, int sizeof_OldScalarFields, struct DNA_OldScalarFields *OldScalarFields, DNA_FLOAT val)
{
  for (int iPoint = 0; iPoint < NumericsFD->NPoints; iPoint++)
  {
    for (int j = 0; j < sizeof_OldScalarFields; j++)
    {
      OldScalarFields->o[j].val[iPoint] = val;
    }
  }

  return 0;
}
