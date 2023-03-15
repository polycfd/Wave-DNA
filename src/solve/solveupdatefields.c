#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Functions to update fields that do not require discretization 
---------------------------------------------------------**/

/** Function to update all old fields required for the time integration **/
int SolveUpdateOldFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  for (int iPoint = 0; iPoint < RunOptions->NumericsFD.NPoints; iPoint++)
  {
    Fields->OldPhiField.Old_phi[1][iPoint] = Fields->OldPhiField.Old_phi[0][iPoint];
    Fields->OldPhiField.Old_phi[0][iPoint] = Fields->PhiField.phi[iPoint];

    Fields->OldPhiField.Old_dXI1_qphi[1][iPoint] = Fields->OldPhiField.Old_dXI1_qphi[0][iPoint];
    Fields->OldPhiField.Old_dXI1_qphi[0][iPoint] = Fields->PhiField.dXI1_qphi[iPoint];

    Fields->OldPhiField.Old_dXI1_phi[1][iPoint] = Fields->OldPhiField.Old_dXI1_phi[0][iPoint];
    Fields->OldPhiField.Old_dXI1_phi[0][iPoint] = Fields->PhiField.dXI1_phi[iPoint];

    Fields->OldPhiField.Old_dXI2_phi[1][iPoint] = Fields->OldPhiField.Old_dXI2_phi[0][iPoint];
    Fields->OldPhiField.Old_dXI2_phi[0][iPoint] = Fields->PhiField.dXI2_phi[iPoint];
  }

  return 0;
}

/** The term qphi represents the product of the acoustic perturbation potential phi
  and the velocity q = dXI/dt of the grid points in the computational domain as viewed
  from a moving point in the physical domain (also taking the Jacobain into ccoutn) **/
int SolveCalcqphi(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  for (int iPoint = 0; iPoint < RunOptions->NumericsFD.NPoints; iPoint++)
  {
    Fields->PhiField.qphi[iPoint] = Fields->Grid.q[iPoint] * Fields->PhiField.phi[iPoint];
  }

  return 0;
}