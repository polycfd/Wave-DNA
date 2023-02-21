#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to allocate memory
---------------------------------------------------------**/

int MemoryAllocGrid (int NPoints, struct DNA_Grid *Grid)
{
  Grid->x = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->XI = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->q = malloc(NPoints * sizeof(DNA_FLOAT));
  Grid->dtq = malloc(NPoints * sizeof(DNA_FLOAT));
  
  return 0;
}

int MemoryAllocScalarField (int NPoints, struct DNA_ScalarField *ScalarField)
{
  ScalarField->val = malloc(NPoints * sizeof(DNA_FLOAT));
    
  return 0;
}

int MemoryAllocOldScalarFields(int NPoints, int structSize, struct DNA_OldScalarFields *OldScalarFields)
{
  int iStruct;
  
  for(iStruct=0; iStruct<structSize; iStruct++)
  {
    (OldScalarFields->o[iStruct]).val = malloc(NPoints * sizeof(DNA_FLOAT));
  }
  
  return 0;
}

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
  Probes->SamplePressure = (DNA_FLOAT **)malloc(rows * sizeof(DNA_FLOAT*));
  for(int i = 0; i < rows; i++) 
  {
    Probes->SamplePressure[i] = (DNA_FLOAT *)malloc(cols * sizeof(DNA_FLOAT));
  }

  return 0;
}

