#include "DNA.h"
#include "DNA-functions.h"



int TransformEqn_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties)
{
  int iPoint;
  
  DNA_FLOAT NLflag = RunOptions->LagrangianDensityFlag;
  
  DNA_FLOAT transformed_dXI1_phi;
  DNA_FLOAT transformed_dXI2_phi;
  DNA_FLOAT transformed_MixedDerivative;
  DNA_FLOAT geometryFactor;
  DNA_FLOAT transformed_dt1;
  DNA_FLOAT transformed_dt2;
  DNA_FLOAT q;
  DNA_FLOAT dqdt;
  DNA_FLOAT A1;
  DNA_FLOAT A2;
  DNA_FLOAT AG;
  DNA_FLOAT AL;
  DNA_FLOAT Agrav;
  DNA_FLOAT B1;
  DNA_FLOAT B2;
  DNA_FLOAT BM;
  DNA_FLOAT BG;
  DNA_FLOAT BL;
  DNA_FLOAT BLr;
  DNA_FLOAT N;
  DNA_FLOAT Uterm;
  DNA_FLOAT UM;
  DNA_FLOAT UG;
  DNA_FLOAT UL;
  DNA_FLOAT u;
  DNA_FLOAT dudx;
  DNA_FLOAT DuDt;

  DNA_FLOAT pow2_c0 = DNA_POW2(FluidProperties->c0);
  DNA_FLOAT K = 2.0*(FluidProperties->beta - 1.0*NLflag)/pow2_c0;
  
  DNA_FLOAT invPow1_dt = (DNA_FLOAT)1.0/RunOptions->NumericsFD.dt;
  DNA_FLOAT invPow2_dt = (DNA_FLOAT)1.0/DNA_POW2(RunOptions->NumericsFD.dt);
  
  DNA_FLOAT a02_times_invPow2_dt = RunOptions->NumericsFD.FDCoeffs.adt_order1.dt2[0]*invPow2_dt;
  DNA_FLOAT a01_times_invPow1_dt = RunOptions->NumericsFD.FDCoeffs.adt_order1.dt1[0]*invPow1_dt;

  DNA_FLOAT detJ = Fields->Grid.detJacobi; 
  DNA_FLOAT ddetJdxi = Fields->Grid.dxidetJacobi;
  DNA_FLOAT divq = Fields->Grid.divq;

  DNA_FLOAT BMblend = 0.5;

  
  for(iPoint=0; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  { 
    q = Fields->Grid.q[iPoint];
    dqdt = Fields->Grid.dtq[iPoint];

    u = Fields->BackgroundVelocity.val[iPoint];
    dudx = Fields->GradBackgroundVelocity.val[iPoint];
    DuDt = Fields->dt1material_BackgroundVelocity.val[iPoint];
    geometryFactor = (*RunOptions->GeometricalDecay)(Fields->Grid.x[iPoint]);
    
    transformed_dXI1_phi = detJ*Fields->dXI1_phi.val[iPoint];    
    transformed_dXI2_phi = detJ*
    (
        ddetJdxi*Fields->dXI1_phi.val[iPoint] 
      + detJ*Fields->dXI2_phi.val[iPoint]
    );

    transformed_MixedDerivative = detJ*
    (
        q*Fields->dXI2_phi.val[iPoint]
      + Fields->dt1_dXI1_phi.val[iPoint]
      + Fields->Grid.divq*Fields->dXI1_phi.val[iPoint]
    );

    transformed_dt1 = Fields->dt1_phi.val[iPoint] + q*Fields->dXI1_phi.val[iPoint];
    transformed_dt2 = 
      Fields->dt2_phi.val[iPoint] 
    + 2.0*q*Fields->dt1_dXI1_phi.val[iPoint]
    + (dqdt + q*divq)*Fields->dXI1_phi.val[iPoint]
    + DNA_POW2(q)*Fields->dXI2_phi.val[iPoint];

    UM = 2.0*u*detJ;
    UG = 2.0*u*detJ*Fields->Grid.divq;
    UL = 2.0*u*detJ*q;

    //Uterm = 2.0*u*transformed_MixedDerivative + DuDt*transformed_dXI1_phi;
    Uterm = 
      UM*Fields->dt1_dXI1_phi.val[iPoint] 
    + UG*Fields->dXI1_phi.val[iPoint]
    + DuDt*transformed_dXI1_phi;

    Agrav = (*RunOptions->WaveCalcGravitationalPotential)(RunOptions, FluidProperties, Fields->Grid.x[iPoint]);
    
    AG =
      dudx*transformed_dXI1_phi*NLflag
    + 2.0*transformed_MixedDerivative*NLflag
    + K*u*Uterm
    + Agrav;

    AL = 
      DNA_POW2(u) - pow2_c0
    + 2.0*u*transformed_dXI1_phi*NLflag
    + K*DNA_POW3(u)*transformed_dXI1_phi;

    A1 = K*(Uterm + DNA_POW2(u)*transformed_dXI2_phi);
    A2 = 1.0 + K*u*transformed_dXI1_phi;

    N = -K*transformed_dt1*transformed_dt2;

    B1 = A1 + K*transformed_dt2;
    B2 = A2 + K*transformed_dt1;

    BM = 2.0*B2;

    BG = 
      B1*q
    + B2*(dqdt + q*divq)
    + AG*detJ
    + AL*detJ*ddetJdxi;

    BL = B2*DNA_POW2(q) + AL*DNA_POW2(detJ) + UL;

    BLr = -pow2_c0*geometryFactor*detJ;
    
    Fields->BB.val[iPoint] =
    - B1*Fields->sum_aphi_dt1.val[iPoint]*invPow1_dt
    - B2*Fields->sum_aphi_dt2.val[iPoint]*invPow2_dt
    - BG*Fields->dXI1_phi.val[iPoint]
    - BM*q*Fields->dt1_dXI1_phi.val[iPoint]*BMblend
    - BM*Fields->dt1_dXI1_qphi.val[iPoint]*(1.0 - BMblend)
    + BM*divq*Fields->sum_aphi_dt1.val[iPoint]*invPow1_dt*(1.0 - BMblend)
    + BM*dqdt*Fields->dXI1_phi.val[iPoint]*(1.0 - BMblend)
    - Uterm
    - N;
    
    Fields->AA.val[iPoint] =
      B1*a01_times_invPow1_dt
    + B2*a02_times_invPow2_dt
    - BM*divq*a01_times_invPow1_dt*(1.0 - BMblend)
    - BM*(2.0*DNA_POW2(divq) + detJ*MovingBoundary->Udot)*(1.0 - BMblend);

    Fields->RHS.val[iPoint] = -BL*Fields->dXI2_phi.val[iPoint] - BLr*Fields->dXI1_phi.val[iPoint];
  }
  
  return 0;
}











int TransformEqn_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties)
{
  int iPoint;
  
  DNA_FLOAT NLflag = RunOptions->LagrangianDensityFlag;
  
  DNA_FLOAT geometryFactor;
  DNA_FLOAT transformed_dXI1_phi;
  DNA_FLOAT transformed_dt1;
  DNA_FLOAT q;
  DNA_FLOAT A2;
  DNA_FLOAT B2;
  DNA_FLOAT AL;
  DNA_FLOAT BL;
  DNA_FLOAT BLr;
  DNA_FLOAT UL;
  DNA_FLOAT u;

  DNA_FLOAT pow2_c0 = DNA_POW2(FluidProperties->c0);
  DNA_FLOAT K = 2.0*(FluidProperties->beta - 1.0*NLflag)/pow2_c0;
  
  DNA_FLOAT detJ = Fields->Grid.detJacobi; 
  
  for(iPoint=0; iPoint<(RunOptions->NumericsFD.NPoints); iPoint++)
  {         
    q = Fields->Grid.q[iPoint];

    geometryFactor = (*RunOptions->GeometricalDecay)(Fields->Grid.x[iPoint]);

    u = Fields->BackgroundVelocity.val[iPoint];

    //Should the old dXI1_phi field be taken?
    transformed_dXI1_phi = detJ*Fields->Old_dXI1_phi.o[0].val[iPoint]; 
    transformed_dt1 = Fields->dt1_phi.val[iPoint] + q*Fields->Old_dXI1_phi.o[0].val[iPoint];

    UL = 2.0*u*detJ*q;

    A2 = 1.0 + K*u*transformed_dXI1_phi;
    
    B2 = A2 + K*transformed_dt1;

    AL = 
      DNA_POW2(u) - pow2_c0
    + 2.0*u*transformed_dXI1_phi*NLflag
    + K*DNA_POW3(u)*transformed_dXI1_phi;

    BL = B2*DNA_POW2(q) + AL*DNA_POW2(detJ) + UL;

    BLr = -pow2_c0*geometryFactor*detJ;
    
    Fields->RHS.val[iPoint] = -BL*Fields->dXI2_phi.val[iPoint] - BLr*Fields->dXI1_phi.val[iPoint];
  }
  
  return 0;
}
