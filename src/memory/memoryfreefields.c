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

  free((Fields->BackgroundFlowField).BackgroundVelocity);
  free((Fields->BackgroundFlowField).GradBackgroundVelocity);
  free((Fields->BackgroundFlowField).dt1_BackgroundVelocity);
  free((Fields->BackgroundFlowField).dt1material_BackgroundVelocity);

  free((Fields->PhiField).phi);
  free((Fields->PhiField).qphi);
  free((Fields->PhiField).dXI2_phi);
  free((Fields->PhiField).dXI1_phi);
  free((Fields->PhiField).dXI1_qphi);
  free((Fields->PhiField).dt1_phi);
  free((Fields->PhiField).dt2_phi);
  free((Fields->PhiField).dt1_dXI1_qphi);
  free((Fields->PhiField).dt1_dXI1_phi);

  free((Fields->PhiField).sum_aphi_dt1);
  free((Fields->PhiField).sum_aphi_dt2);
  free((Fields->PhiField).sum_adXI1_phi_dt1);
  free((Fields->PhiField).sum_adXI1_qphi_dt1);

  free((Fields->PhiField).RHS);
  free((Fields->PhiField).AA);
  free((Fields->PhiField).BB);
  free((Fields->PhiField).BBbyAA);
  free((Fields->PhiField).phi1_initGuess);

  free((Fields->PhiField).PressureField);

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
