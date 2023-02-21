#include "DNA.h"
#include "DNA-functions.h"


int SolveCorrectLaplacian(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDGradient(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI1_phi);
  FDLaplacian(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI2_phi);
  
  return 0;
}



int SolveExplicit_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;
  
  for(iPoint=1; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->phi.val[iPoint] =
    (Fields->RHS.val[iPoint] + Fields->BB.val[iPoint])/Fields->AA.val[iPoint];
  }
  
  return 0;
}



int SolveExplicit_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;
  DNA_FLOAT gamma = RunOptions->NumericsFD.correctorWeight;
  
  for(iPoint=1; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {    
    Fields->phi.val[iPoint] =
      (1.0 - gamma)*Fields->phi1_initGuess.val[iPoint]
    + gamma*Fields->BBbyAA.val[iPoint]
    + gamma*Fields->RHS.val[iPoint]/Fields->AA.val[iPoint]; 
  }
  
  return 0;
}



int SolveStoreInitialGuess(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;
  
  for(iPoint=0; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BBbyAA.val[iPoint] = Fields->BB.val[iPoint]/(Fields->AA.val[iPoint]);
    Fields->phi1_initGuess.val[iPoint] = Fields->phi.val[iPoint]; 
  }
  
  return 0;
}


