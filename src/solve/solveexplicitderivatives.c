#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions explicitly compute the required derivatives in the coordinates 
(XI,t) of the computational domain by calling the functions in fddx.c and fddt.c
---------------------------------------------------------**/

int SolveExplicitDerivatives_dx(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDGradient(&RunOptions->NumericsFD, Fields->PhiField.phi, Fields->PhiField.dXI1_phi);

  FDLaplacian(&RunOptions->NumericsFD, Fields->PhiField.phi, Fields->PhiField.dXI2_phi);

  FDGradient(&RunOptions->NumericsFD, Fields->PhiField.qphi, Fields->PhiField.dXI1_qphi);

  return 0;
}

int SolveExplicitDerivatives_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  // First partial time derivative of phi
  FDdt1(&RunOptions->NumericsFD, Fields->PhiField.sum_aphi_dt1, Fields->PhiField.phi, Fields->PhiField.dt1_phi);

  // Second partial time derivative of phi
  FDdt2(&RunOptions->NumericsFD, Fields->PhiField.sum_aphi_dt2, Fields->PhiField.phi, Fields->PhiField.dt2_phi);

  // First partial time derivative of dXI1_phi (in coordinates of the computational domain)
  FDdt1(&RunOptions->NumericsFD, Fields->PhiField.sum_adXI1_phi_dt1, Fields->PhiField.dXI1_phi, Fields->PhiField.dt1_dXI1_phi);

  // First partial time derivative of dt1_dXI1_qphi (in coordinates of the computational domain)
  FDdt1(&RunOptions->NumericsFD, Fields->PhiField.sum_adXI1_qphi_dt1, Fields->PhiField.dXI1_qphi, Fields->PhiField.dt1_dXI1_qphi);

  return 0;
}
