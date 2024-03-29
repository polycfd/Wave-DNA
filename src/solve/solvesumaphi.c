#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
For a finite difference scheme for the k-th partial time derivative
d^k phi/dt^k = a0*phi_new + Sum((a*phi)_old)/dt^k, the following function calls the
functions to compute Sum((a*phi)_old).
---------------------------------------------------------**/

int SolveSumFDcoeffTimesPhi_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDSumFDcoeffTimesPhi_dt1(&RunOptions->NumericsFD, Fields->PhiField.sum_aphi_dt1, Fields->OldPhiField.Old_phi);
  FDSumFDcoeffTimesPhi_dt2(&RunOptions->NumericsFD, Fields->PhiField.sum_aphi_dt2, Fields->OldPhiField.Old_phi);
  FDSumFDcoeffTimesPhi_dt1(&RunOptions->NumericsFD, Fields->PhiField.sum_adXI1_phi_dt1, Fields->OldPhiField.Old_dXI1_phi);
  FDSumFDcoeffTimesPhi_dt1(&RunOptions->NumericsFD, Fields->PhiField.sum_adXI1_qphi_dt1, Fields->OldPhiField.Old_dXI1_qphi);

  return 0;
}
