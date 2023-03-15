#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
This function calls all the routines for memory allocation as implemented in
<memoryalloc.c> and it is call upon initialization (see <initializesimulation.c>). 
---------------------------------------------------------**/

int MemoryAllocFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  if (RunOptions->Probes.nSamplePoints > 0)
  {
    MemoryAllocSampleVectors(RunOptions->writeFrequency, &RunOptions->Probes);
  }

  MemoryAllocGrid(RunOptions->NumericsFD.NPoints, &Fields->Grid);
  MemoryAllocBackgroundFlowField(RunOptions->NumericsFD.NPoints, &Fields->BackgroundFlowField);
  MemoryAllocPhiField(RunOptions->NumericsFD.NPoints, &Fields->PhiField);

  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->phi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->qphi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dXI2_phi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dXI1_phi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dXI1_qphi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dt1_phi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dt2_phi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dt1_dXI1_qphi);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->dt1_dXI1_phi);

  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->sum_aphi_dt1);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->sum_aphi_dt2);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->sum_adXI1_phi_dt1);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->sum_adXI1_qphi_dt1);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->phi1_initGuess);

  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->RHS);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->AA);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->BB);
  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->BBbyAA);

  MemoryAllocScalarField(RunOptions->NumericsFD.NPoints, &Fields->PressureField);

  MemoryAllocOldScalarFields(RunOptions->NumericsFD.NPoints, Fields->sizeof_OldScalarFields, &Fields->Old_phi);
  MemoryAllocOldScalarFields(RunOptions->NumericsFD.NPoints, Fields->sizeof_OldScalarFields, &Fields->Old_dXI1_phi);
  MemoryAllocOldScalarFields(RunOptions->NumericsFD.NPoints, Fields->sizeof_OldScalarFields, &Fields->Old_dXI2_phi);
  MemoryAllocOldScalarFields(RunOptions->NumericsFD.NPoints, Fields->sizeof_OldScalarFields, &Fields->Old_dXI1_qphi);

  return 0;
}
