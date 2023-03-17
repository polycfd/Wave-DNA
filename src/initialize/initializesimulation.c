#include <sys/stat.h>

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
  // Generate results directory if it does not yet exist.
  struct stat st = {0};
  if (stat("./results", &st) == -1) mkdir("results", 0700);

  RunOptions->writeCount = 0;
  IOWriteHeader(RunOptions);

  RunOptions->pExcitation = 0;

  MovingBoundary->R = MovingBoundary->R0;
  MovingBoundary->U = 0.0;
  MovingBoundary->Udot = 0.0;

  FDFiniteDifferenceCoeffs(&RunOptions->NumericsFD);

  MemoryAllocFields(RunOptions, Fields);

  GridComputationalDomain(&RunOptions->NumericsFD, &Fields->Grid);
  Fields->Grid.xmov = MovingBoundary->R;
  GridPhysicalDomain(&RunOptions->NumericsFD, &Fields->Grid);
  GridMotion(RunOptions, Fields, MovingBoundary);

  for (int iPoint = 0; iPoint < RunOptions->NumericsFD.NPoints; iPoint++)
  {
    Fields->Grid.q[iPoint] = 0.0;
    Fields->Grid.dtq[iPoint] = 0.0;

    Fields->BackgroundFlowField.BackgroundVelocity[iPoint] = 0.0;
    Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] = 0.0;
    Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] = 0.0;
    Fields->BackgroundFlowField.dt1material_BackgroundVelocity[iPoint] = 0.0;
  }

  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.phi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.qphi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dXI1_phi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dXI1_qphi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dXI2_phi, 0.0);

  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dt1_phi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dt2_phi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dt1_dXI1_phi, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.dt1_dXI1_qphi, 0.0);

  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.PressureField, 0.0);

  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.phi1_initGuess, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.RHS, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.AA, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.BB, 0.0);
  InitializeScalarField(&RunOptions->NumericsFD, Fields->PhiField.BBbyAA, 0.0);

  InitializeOldPhiFields(&RunOptions->NumericsFD, Fields->OldPhiField.Old_phi, 0.0);
  InitializeOldPhiFields(&RunOptions->NumericsFD, Fields->OldPhiField.Old_dXI1_phi, 0.0);
  InitializeOldPhiFields(&RunOptions->NumericsFD, Fields->OldPhiField.Old_dXI2_phi, 0.0);
  InitializeOldPhiFields(&RunOptions->NumericsFD, Fields->OldPhiField.Old_dXI1_qphi, 0.0);

  SolveSumFDcoeffTimesPhi_dt(RunOptions, Fields);

  IOIdentifySampleIDs(RunOptions, Fields);

  return 0;
}
