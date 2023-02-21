#include "DNA.h"
#include "DNA-functions.h"


int SolveExplicitDerivatives_dx(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{   
  // dXI1_phi
  FDGradient(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI1_phi);

  // dXI2_phi
  FDLaplacian(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI2_phi);

  FDGradient(&RunOptions->NumericsFD, &Fields->qphi, &Fields->dXI1_qphi);

  return 0;
}




int SolveExplicitDerivatives_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{ 
  // First time derivative of phi1
  FDdt1(&RunOptions->NumericsFD, &Fields->sum_aphi_dt1, &Fields->phi, &Fields->dt1_phi);

  // Second time derivative of phi1
  FDdt2(&RunOptions->NumericsFD, &Fields->sum_aphi_dt2, &Fields->phi, &Fields->dt2_phi);

  // First time derivative of dXI1_phi
  FDdt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_phi_dt1, &Fields->dXI1_phi, &Fields->dt1_dXI1_phi);

  FDdt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_qphi_dt1, &Fields->dXI1_qphi, &Fields->dt1_dXI1_qphi);

  return 0;
}
