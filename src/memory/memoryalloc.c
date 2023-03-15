#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to allocate memory 
---------------------------------------------------------**/

int MemoryAllocGrid(int NPoints, struct DNA_Grid *Grid)
{
  Grid->x = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->XI = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->q = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->dtq = malloc(NPoints * sizeof(DNA_FLOAT));

  return 0;
}

int MemoryAllocBackgroundFlowField(int NPoints, struct DNA_BackgroundFlowField *BackgroundFlowField)
{
  BackgroundFlowField->BackgroundVelocity = malloc(NPoints * sizeof(DNA_FLOAT));
  BackgroundFlowField->GradBackgroundVelocity = malloc(NPoints * sizeof(DNA_FLOAT));
  BackgroundFlowField->dt1_BackgroundVelocity = malloc(NPoints * sizeof(DNA_FLOAT));
  BackgroundFlowField->dt1material_BackgroundVelocity = malloc(NPoints * sizeof(DNA_FLOAT));

  return 0;
}

int MemoryAllocPhiField(int NPoints, struct DNA_PhiField *PhiField)
{
  PhiField->phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->qphi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dXI1_phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dXI2_phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dXI1_qphi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dt1_phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dt2_phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dt1_dXI1_phi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->dt1_dXI1_qphi = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->PressureField = malloc(NPoints * sizeof(DNA_FLOAT));

  PhiField->sum_aphi_dt1 = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->sum_aphi_dt2 = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->sum_adXI1_phi_dt1 = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->sum_adXI1_qphi_dt1 = malloc(NPoints * sizeof(DNA_FLOAT));

  PhiField->phi1_initGuess = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->RHS = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->AA = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->BB = malloc(NPoints * sizeof(DNA_FLOAT));
  PhiField->BBbyAA = malloc(NPoints * sizeof(DNA_FLOAT));

  return 0;
}

int MemoryAllocOldPhiField(int NPoints, struct DNA_OldPhiField *OldPhiField)
{
  OldPhiField->Old_phi = (DNA_FLOAT **) malloc((int) 2 * sizeof(DNA_FLOAT *));
  OldPhiField->Old_phi[0] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));
  OldPhiField->Old_phi[1] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));

  OldPhiField->Old_dXI1_phi = (DNA_FLOAT **) malloc((int) 2 * sizeof(DNA_FLOAT *));
  OldPhiField->Old_dXI1_phi[0] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));
  OldPhiField->Old_dXI1_phi[1] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));

  OldPhiField->Old_dXI2_phi = (DNA_FLOAT **) malloc((int) 2 * sizeof(DNA_FLOAT *));
  OldPhiField->Old_dXI2_phi[0] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));
  OldPhiField->Old_dXI2_phi[1] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));

  OldPhiField->Old_dXI1_qphi = (DNA_FLOAT **) malloc((int) 2 * sizeof(DNA_FLOAT *));
  OldPhiField->Old_dXI1_qphi[0] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));
  OldPhiField->Old_dXI1_qphi[1] = (DNA_FLOAT *) malloc(NPoints * sizeof(DNA_FLOAT));

  return 0;
}

/**--------------------------------------------------------- 
The following two functions allocate memory for the acoustic pressure probes,
taken at the grid point IDs that correspond to the sample point locations
specified in the options file.

The function <MemoryAllocSamplePoints> allocates sample IDs (int) and the
sample points (float). The function is called in the <ioreadoptions.c> file,
where NPoints is the number of sample points as specificed in the options file.

Subsequently, <MemoryAllocSampleVectors> is called in <MemoryAllocFields>
(see <memoryallocfields.c>), which again is called in <InitializeSimulation>
(see <initializesimulation.c>) in order to allocate the memory per sample point
based on the write frequency as specified in the options file (<writeFrequency>).
The write frequency is an integer indicating the number of subsequent time steps
after which the field data and, on this occasion, the sample data is written to
disk. 
---------------------------------------------------------**/

int MemoryAllocSamplePoints(int NPoints, struct DNA_Probes *Probes)
{
  Probes->SamplePoints = malloc(NPoints * sizeof(DNA_FLOAT));
  Probes->SampleIDs = malloc(NPoints * sizeof(int));

  return 0;
}

int MemoryAllocSampleVectors(int writeFrequency, struct DNA_Probes *Probes)
{
  Probes->SampleTime = malloc(writeFrequency * sizeof(DNA_FLOAT));
  Probes->TimeIDs = malloc(writeFrequency * sizeof(int));

  int rows = writeFrequency;
  int cols = Probes->nSamplePoints;
  Probes->SamplePressure = (DNA_FLOAT **) malloc(rows * sizeof(DNA_FLOAT *));

  for (int i = 0; i < rows; i++)
  {
    Probes->SamplePressure[i] = (DNA_FLOAT *) malloc(cols * sizeof(DNA_FLOAT));
  }

  return 0;
}
