#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
This function frees all alocated memory upon termination of the simulation
and it is called in main.c.
---------------------------------------------------------**/

int MemoryFreeFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  if (RunOptions->Probes.nSamplePoints > 0)
  {
    free(RunOptions->Probes.SamplePoints);
    MemoryFreeSampleVectors(RunOptions->writeFrequency, &RunOptions->Probes);
    free(RunOptions->Probes.SampleIDs);
  }

  free((Fields->Grid).x);
  free((Fields->Grid).XI);
  free((Fields->Grid).q);
  free((Fields->Grid).dtq);

  free((Fields->phi).val);
  free((Fields->qphi).val);
  free((Fields->dXI2_phi).val);
  free((Fields->dXI1_phi).val);
  free((Fields->dXI1_qphi).val);
  free((Fields->dt1_phi).val);
  free((Fields->dt2_phi).val);
  free((Fields->dt1_dXI1_qphi).val);
  free((Fields->dt1_dXI1_phi).val);

  free((Fields->sum_aphi_dt1).val);
  free((Fields->sum_aphi_dt2).val);
  free((Fields->sum_adXI1_phi_dt1).val);
  free((Fields->sum_adXI1_qphi_dt1).val);

  free((Fields->RHS).val);
  free((Fields->AA).val);
  free((Fields->BB).val);
  free((Fields->BBbyAA).val);
  free((Fields->phi1_initGuess).val);

  free((Fields->BackgroundVelocity).val);
  free((Fields->GradBackgroundVelocity).val);
  free((Fields->dt1_BackgroundVelocity).val);
  free((Fields->dt1material_BackgroundVelocity).val);

  free((Fields->PressureField).val);

  MemoryFreeOldScalarFields(Fields->sizeof_OldScalarFields, &Fields->Old_phi);
  MemoryFreeOldScalarFields(Fields->sizeof_OldScalarFields, &Fields->Old_dXI1_phi);
  MemoryFreeOldScalarFields(Fields->sizeof_OldScalarFields, &Fields->Old_dXI2_phi);
  MemoryFreeOldScalarFields(Fields->sizeof_OldScalarFields, &Fields->Old_dXI1_qphi);

  return 0;
}

int MemoryFreeOldScalarFields(int structSize, struct DNA_OldScalarFields *OldScalarFields)
{
  for (int iStruct = 0; iStruct < structSize; iStruct++)
  {
    free((OldScalarFields->o[iStruct]).val);
  }

  return 0;
}

int MemoryFreeSampleVectors(int writeFrequency, struct DNA_Probes *Probes)
{
  for (int i = 0; i < writeFrequency; i++)
  {
    free(Probes->SamplePressure[i]);
  }

  free(Probes->SampleTime);
  free(Probes->TimeIDs);
  free(Probes->SamplePressure);

  return 0;
}
