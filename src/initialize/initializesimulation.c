#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
This function calls the functions to allocate memory for the fields and
assigns initial values to all variables and fields. Furthermore, the
function IOWriteHeader is called to create the output files and to write
the file headers. If an output file already exists, the simulation data 
will be appended to the existing data content. 

Furthermore, the sample IDs are identified based on the locations
specified by the user. The acoustic pressure at those IDs is written to
the output file probes.dat. Note that the corresponding grid points may
be moving.
---------------------------------------------------------**/

int InitializeSimulation(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{
  int iPoint;
  
  RunOptions->writeCount = 0;
  IOWriteHeader(RunOptions);
  
  RunOptions->pExcitation = 0;

  MovingBoundary->R = MovingBoundary->R0;
  MovingBoundary->U = 0.0;
  MovingBoundary->Udot = 0.0;
  
  FDFiniteDifferenceCoeffs(&RunOptions->NumericsFD); 

  Fields->sizeof_OldScalarFields = sizeof(struct DNA_OldScalarFields)/sizeof(DNA_FLOAT);
  
  MemoryAllocFields(RunOptions, Fields); 
  
  GridComputationalDomain(&(RunOptions->NumericsFD), &(Fields->Grid));
  Fields->Grid.xmov = MovingBoundary->R;
  GridPhysicalDomain(&(RunOptions->NumericsFD), &(Fields->Grid));
  GridMotion(RunOptions, Fields, MovingBoundary);
  
  for (iPoint=0; iPoint<RunOptions->NumericsFD.NPoints; iPoint++)
  {
    Fields->Grid.q[iPoint] = 0.0;
    Fields->Grid.dtq[iPoint] = 0.0;
  }
  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->phi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->qphi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dXI1_phi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dXI1_qphi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dXI2_phi), 0.0);
  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt1_phi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt2_phi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt1_dXI1_phi), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt1_dXI1_qphi), 0.0);

  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->BackgroundVelocity), 0.0);  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->GradBackgroundVelocity), 0.0);  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt1_BackgroundVelocity), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->dt1material_BackgroundVelocity), 0.0);
  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->PressureField), 0.0);
  
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->phi1_initGuess), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->RHS), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->AA), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->BB), 0.0);
  InitializeConstScalarField(&(RunOptions->NumericsFD), &(Fields->BBbyAA), 0.0);
  
  InitializeConstOldScalarFields(&(RunOptions->NumericsFD), Fields->sizeof_OldScalarFields, &(Fields->Old_phi), 0.0);
  InitializeConstOldScalarFields(&(RunOptions->NumericsFD), Fields->sizeof_OldScalarFields, &(Fields->Old_dXI1_qphi), 0.0);
  InitializeConstOldScalarFields(&(RunOptions->NumericsFD), Fields->sizeof_OldScalarFields, &(Fields->Old_dXI1_phi), 0.0);
  InitializeConstOldScalarFields(&(RunOptions->NumericsFD), Fields->sizeof_OldScalarFields, &(Fields->Old_dXI2_phi), 0.0);

  SolveSumAy_dt(RunOptions, Fields);
  
  IOIdentifySampleIDs(RunOptions, Fields);
  
  return 0;
}
