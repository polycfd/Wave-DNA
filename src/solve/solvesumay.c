#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
For a finite difference scheme for the k-th partial time derivative
d^k phi/dt^k = a0*phi_new + Sum((a*phi)_old)/dt^k, the following function calls the
functions to compute Sum((a*phi)_old).
---------------------------------------------------------**/

int SolveSumAy_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDSumAY_dt1(&RunOptions->NumericsFD, &Fields->sum_aphi_dt1, &Fields->Old_phi);
  FDSumAY_dt2(&RunOptions->NumericsFD, &Fields->sum_aphi_dt2, &Fields->Old_phi);
  FDSumAY_dt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_phi_dt1, &Fields->Old_dXI1_phi);
  FDSumAY_dt1(&RunOptions->NumericsFD, &Fields->sum_adXI1_qphi_dt1, &Fields->Old_dXI1_qphi);

  return 0;
}
