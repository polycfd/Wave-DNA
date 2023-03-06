#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following two functions convert the excitation pressure to the excitation
potential at the excitation node <iExcitation>
(<TransformExcitationPressureToExcitationPotential>) and the acoustic perturbation
potential field phi into the acoustic pressure field
(<TransformPotentialFieldToPressureField>). The variable names follow a similar
pattern as in <transformeqn.c>.

The conversion essentially relies on the transformed form of the transient
Bernoulli equation (expressed in terms of phi). After the terms related to phi are
constructed, the acoustic pressure follows by forward substitution, whereas the
excitation potential is the "new" term in the time-marching scheme (corresponding
to the a0 coefficient of the time integration scheme).
---------------------------------------------------------**/

int TransformExcitationPressureToExcitationPotential(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties)
{
  DNA_FLOAT NLflag = 0.0;//RunOptions->LagrangianDensityFlag;
  
  DNA_FLOAT phi;
  
  DNA_FLOAT p = RunOptions->pExcitation;
  DNA_FLOAT b = Fields->sum_aphi_dt1.val[RunOptions->iExcitation];
  DNA_FLOAT dt = RunOptions->NumericsFD.dt;
  DNA_FLOAT a0 = RunOptions->NumericsFD.FDCoeffs.adt_order1.dt1[0];
  DNA_FLOAT detJ = Fields->Grid.detJacobi;
  DNA_FLOAT dXI1_phi = Fields->dXI1_phi.val[RunOptions->iExcitation];
  DNA_FLOAT u0 = Fields->BackgroundVelocity.val[RunOptions->iExcitation];
  DNA_FLOAT q = Fields->Grid.q[RunOptions->iExcitation];
  DNA_FLOAT dtphi = Fields->dt1_phi.val[RunOptions->iExcitation];
  
  DNA_FLOAT transformed_dXI1_phi = detJ*dXI1_phi;
  DNA_FLOAT transformed_dt1 = dtphi + q*dXI1_phi;
  
  DNA_FLOAT rho0 = FluidProperties->rho0;
  DNA_FLOAT pow2_c0 = DNA_POW2(FluidProperties->c0);
    
  DNA_FLOAT A1 = 1.0 - (transformed_dt1 + u0*transformed_dXI1_phi)/pow2_c0*NLflag;
  DNA_FLOAT AG = u0 + 0.5*(1.0 - DNA_POW2(u0)/pow2_c0)*transformed_dXI1_phi*NLflag;
  DNA_FLOAT N = (DNA_FLOAT)0.5*DNA_POW2(transformed_dt1)/pow2_c0*NLflag;
  
  phi = -b/a0 -dt/(a0*A1)*((A1*q + AG*detJ)*dXI1_phi + N + p/rho0);
  
  Fields->phi.val[RunOptions->iExcitation] = phi;  
  
  return 0;
}

int TransformPotentialFieldToPressureField(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties)
{
  DNA_FLOAT NLflag = RunOptions->LagrangianDensityFlag;
  
  int iPoint;
  
  DNA_FLOAT p;
  DNA_FLOAT detJ = Fields->Grid.detJacobi;
  
  DNA_FLOAT rho0 = FluidProperties->rho0;
  DNA_FLOAT pow2_c0 = DNA_POW2(FluidProperties->c0);
  
  for (iPoint=1; iPoint<RunOptions->NumericsFD.NPoints; iPoint++)
  {
    DNA_FLOAT dXI1_phi = Fields->dXI1_phi.val[iPoint];
    DNA_FLOAT u0 = Fields->BackgroundVelocity.val[iPoint];
    DNA_FLOAT q = Fields->Grid.q[iPoint];
    DNA_FLOAT dtphi = Fields->dt1_phi.val[iPoint];
  
    DNA_FLOAT transformed_dXI1_phi = detJ*dXI1_phi;
    DNA_FLOAT transformed_dt1 = dtphi + q*dXI1_phi;

    DNA_FLOAT A1 = 1.0 - (transformed_dt1 + u0*transformed_dXI1_phi)/pow2_c0*NLflag;
    DNA_FLOAT AG = u0 + 0.5*(1.0 - DNA_POW2(u0)/pow2_c0)*transformed_dXI1_phi*NLflag;
    DNA_FLOAT N = (DNA_FLOAT)0.5*DNA_POW2(transformed_dt1)/pow2_c0*NLflag;
    
    p = -rho0*(A1*dtphi + (A1*q + AG*detJ)*dXI1_phi + N);
    
    Fields->PressureField.val[iPoint] = p;
  }
  
  Fields->PressureField.val[RunOptions->iExcitation] = RunOptions->pExcitation;
  
  return 0;
}