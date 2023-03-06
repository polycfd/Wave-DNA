#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions explicitly compute the required derivatives in the coordinates 
(XI,t) of the computational domain by calling the functions in fddx.c and fddt.c
---------------------------------------------------------**/

int SolveExplicitDerivatives_dx(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDGradient(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI1_phi);

  FDLaplacian(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI2_phi);

  FDGradient(&RunOptions->NumericsFD, &Fields->qphi, &Fields->dXI1_qphi);

  return 0;
}

int SolveExplicitDerivatives_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  // First partial time derivative of phi
  FDdt1(&RunOptions->NumericsFD, &Fields->sum_aphi_dt1, &Fields->phi, &Fields->dt1_phi);

  // Second partial time derivative of phi
  FDdt2(&RunOptions->NumericsFD, &Fields->sum_aphi_dt2, &Fields->phi, &Fields->dt2_phi);

  // First partial time derivative of dXI1_phi (in coordinates of the computational domain)
  FDdt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_phi_dt1, &Fields->dXI1_phi, &Fields->dt1_dXI1_phi);

  // First partial time derivative of dt1_dXI1_qphi (in coordinates of the computational domain)
  FDdt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_qphi_dt1, &Fields->dXI1_qphi, &Fields->dt1_dXI1_qphi);

  return 0;
}
