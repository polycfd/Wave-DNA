#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to update fields that do not require discretization 
---------------------------------------------------------**/

/** Function to update all old fields required for the time integration **/
int SolveUpdateOldFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;

  for (iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
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

/** The term qphi represents the product of the acoustic perturbation potential phi
  and the velocity q = dXI/dt of the grid points in the computational domain as viewed
  from a moving point in the physical domain (also taking the Jacobain into ccoutn) **/
int SolveCalcqphi(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;

  for (iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->qphi.val[iPoint] = Fields->Grid.q[iPoint] * Fields->phi.val[iPoint];
  }

  return 0;
}