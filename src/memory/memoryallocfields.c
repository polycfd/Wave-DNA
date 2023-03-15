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
  MemoryAllocOldPhiField(RunOptions->NumericsFD.NPoints, &Fields->OldPhiField);

  return 0;
}
