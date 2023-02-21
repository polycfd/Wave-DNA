#include "DNA.h"
#include "DNA-functions.h"

int SolveUpdateOldFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;
  
  for(iPoint=0; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {
    (Fields->Old_phi).o[1].val[iPoint] = (Fields->Old_phi).o[0].val[iPoint];
    (Fields->Old_phi).o[0].val[iPoint] = (Fields->phi).val[iPoint];
    
    (Fields->Old_dXI1_qphi).o[1].val[iPoint] = (Fields->Old_dXI1_qphi).o[0].val[iPoint];
    (Fields->Old_dXI1_qphi).o[0].val[iPoint] = (Fields->dXI1_qphi).val[iPoint];

    (Fields->Old_dXI1_phi).o[1].val[iPoint] = (Fields->Old_dXI1_phi).o[0].val[iPoint];
    (Fields->Old_dXI1_phi).o[0].val[iPoint] = (Fields->dXI1_phi).val[iPoint];
    
    (Fields->Old_dXI2_phi).o[1].val[iPoint] = (Fields->Old_dXI2_phi).o[0].val[iPoint];
    (Fields->Old_dXI2_phi).o[0].val[iPoint] = (Fields->dXI2_phi).val[iPoint];
  }
  
  return 0;
}


int SolveCalcqphi(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;

  for(iPoint=0; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->qphi.val[iPoint] = Fields->Grid.q[iPoint]*Fields->phi.val[iPoint];
  }  

  return 0;
}