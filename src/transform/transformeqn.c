#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following two functions are to transform the governing wave equation for the
predictor step (<TransformEqn_Predictor>) and the corrector step 
(<TransformEqn_Corrector>). 

The transformation into the computational coordinates (XI,t) is done in two steps:

1.) Expansion of the material derivatives in the physical coordinates (x,t) so that
the governing equation is written in terms of partial derivatives only. This yields
a set of coefficients A...

2.) Change of variables from physical coordinates (x,t) into comuptational coordinates
(XI,t): This step generates new terms (chain rule) and the equation obtained in this way
is written in terms of a set of coefficients B... which again depend on the coefficients
A...

The second step further involves a Newton linearization of the nonlinear terms as well
as a special treatment of the mixed derivative, which is written as a blending of
two analytically equivalent terms, but which bahave, when discretized, in such a way
that the numerical solution is not biased towards dispersion or diffusion in the
presence of a moving domain boundary. 

The original form of the mixed term is

2*B2*q*[partial^2 phi/partial xi partial t]

and the equivalent expression is obtained by expanding the expression

(2*B2*)*[partial^2(q*phi))/partial xi partial t]

and solving for the original form. The original form and the equivalent form are weighted
by the blending factors BMblend and (1 - BMblend), repectively. The optimal value of 
BMblend is found to be 0.5 and therefore hard-coded, bit nevertheless, BMblend is
formally included as a local variable in order to clearly identify the terms originating
from this procedure.  
---------------------------------------------------------**/

int TransformEqn_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                           struct DNA_FluidProperties *FluidProperties)
{
  /** This flag is equal to 0 if local nonlinearities are not taken into account
  (optins file setting <LocalNonlinearity> set to <False>) and 1 otherwise 
  (optins file setting <LocalNonlinearity> set to <True>). **/
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

  /** The following variables are not space-dependent. **/

  DNA_FLOAT pow2_c0 = DNA_POW2(FluidProperties->c0);
  DNA_FLOAT K = 2.0 * (FluidProperties->beta - 1.0 * NLflag) / pow2_c0;

  DNA_FLOAT invPow1_dt = (DNA_FLOAT) 1.0 / RunOptions->NumericsFD.dt;
  DNA_FLOAT invPow2_dt = (DNA_FLOAT) 1.0 / DNA_POW2(RunOptions->NumericsFD.dt);

  DNA_FLOAT a02_times_invPow2_dt = RunOptions->NumericsFD.FDCoeffs.adt_order1.dt2[0] * invPow2_dt;
  DNA_FLOAT a01_times_invPow1_dt = RunOptions->NumericsFD.FDCoeffs.adt_order1.dt1[0] * invPow1_dt;

  /** The Jacobian determinan detJ and derivatives thereof are time-dependent if
  the physical domain has a moving boundary **/
  DNA_FLOAT detJ = Fields->Grid.detJacobi;
  DNA_FLOAT ddetJdxi = Fields->Grid.dxidetJacobi;
  /** spatial derivative w.r.t. to XI of the velocity q = dXI/dt of the grid points
  in the computational domain as viewed from a moving point in the physical domain 
  (also taking the Jacobain into ccoutn). Due to the uniform spatial distribution of
  the grid points, this quantity is not space-dependent, even though q is
  space-dependent. **/
  DNA_FLOAT divq = Fields->Grid.divq;

  /** Blending facotor explained above **/
  DNA_FLOAT BMblend = 0.5;

  /** In the coefficient names, "1" and "2" refers to the order of the partial
  time-derivatives to which the coefficients belong, and "G", "L", and "M" means
  "gradient", "Laplacian", and "Mixed" (derivative), respectively. **/
  for (int iPoint = 0; iPoint < RunOptions->NumericsFD.NPoints; iPoint++)
  {
    q = Fields->Grid.q[iPoint];
    dqdt = Fields->Grid.dtq[iPoint];

    /** The derivatives of the backgroiund flow velocity field are expressed
    in terms of the physical coordinates (x,t) and not transformed! This is convenient
    as these expressions are known as analytical expressions in the present modeling
    framework so that there is no neeed to discretize them. Note that the coordinate
    transformation between physical and computational coordinates conveys a one-to-one
    mapping, so that this approach is formally valid. The same applied to the geometry
    factor that accounts for the surface area of the curved wavefront due to the 
    flow geometry. **/
    u = Fields->BackgroundVelocity.val[iPoint];
    dudx = Fields->GradBackgroundVelocity.val[iPoint];
    DuDt = Fields->dt1material_BackgroundVelocity.val[iPoint];
    geometryFactor = RunOptions->GeometricalDecay(Fields->Grid.x[iPoint]);

    /** The prefix "transformed_" indicates that the corresponding terms are
    discretized ones but treated explicitly, either because they only involve spatial
    derivatives, or because they invovle time derivatives, which, however, are treated
    explicitly as part of the coefficients of the transformed wave equation. Therefore,
    there is no need to split the corresponding discrete approximations for the time
    integration scheme so that they can be constructed rightaway. **/
    transformed_dXI1_phi = detJ * Fields->dXI1_phi.val[iPoint];
    transformed_dXI2_phi = detJ * (ddetJdxi * Fields->dXI1_phi.val[iPoint] + detJ * Fields->dXI2_phi.val[iPoint]);

    transformed_MixedDerivative =
        detJ * (q * Fields->dXI2_phi.val[iPoint] + Fields->dt1_dXI1_phi.val[iPoint] + Fields->Grid.divq * Fields->dXI1_phi.val[iPoint]);

    transformed_dt1 = Fields->dt1_phi.val[iPoint] + q * Fields->dXI1_phi.val[iPoint];
    transformed_dt2 = Fields->dt2_phi.val[iPoint] + 2.0 * q * Fields->dt1_dXI1_phi.val[iPoint] + (dqdt + q * divq) * Fields->dXI1_phi.val[iPoint] +
                      DNA_POW2(q) * Fields->dXI2_phi.val[iPoint];

    UM = 2.0 * u * detJ;
    UG = 2.0 * u * detJ * Fields->Grid.divq;
    UL = 2.0 * u * detJ * q;

    Uterm = UM * Fields->dt1_dXI1_phi.val[iPoint] + UG * Fields->dXI1_phi.val[iPoint] + DuDt * transformed_dXI1_phi;

    Agrav = RunOptions->WaveCalcGravitationalPotential(RunOptions, FluidProperties, Fields->Grid.x[iPoint]);

    AG = dudx * transformed_dXI1_phi * NLflag + 2.0 * transformed_MixedDerivative * NLflag + K * u * Uterm + Agrav;

    AL = DNA_POW2(u) - pow2_c0 + 2.0 * u * transformed_dXI1_phi * NLflag + K * DNA_POW3(u) * transformed_dXI1_phi;

    A1 = K * (Uterm + DNA_POW2(u) * transformed_dXI2_phi);
    A2 = 1.0 + K * u * transformed_dXI1_phi;

    N = -K * transformed_dt1 * transformed_dt2;

    B1 = A1 + K * transformed_dt2;
    B2 = A2 + K * transformed_dt1;

    BM = 2.0 * B2;

    BG = B1 * q + B2 * (dqdt + q * divq) + AG * detJ + AL * detJ * ddetJdxi;

    BL = B2 * DNA_POW2(q) + AL * DNA_POW2(detJ) + UL;

    /** Gradient contribution to the Laplacian for curved wavefronts **/
    BLr = -pow2_c0 * geometryFactor * detJ;

    /** Assembling the terms for the predictor step (BB, AA, RHS) **/
    Fields->BB.val[iPoint] =
        -B1 * Fields->sum_aphi_dt1.val[iPoint] * invPow1_dt - B2 * Fields->sum_aphi_dt2.val[iPoint] * invPow2_dt - BG * Fields->dXI1_phi.val[iPoint] -
        BM * q * Fields->dt1_dXI1_phi.val[iPoint] * BMblend - BM * Fields->dt1_dXI1_qphi.val[iPoint] * (1.0 - BMblend) +
        BM * divq * Fields->sum_aphi_dt1.val[iPoint] * invPow1_dt * (1.0 - BMblend) + BM * dqdt * Fields->dXI1_phi.val[iPoint] * (1.0 - BMblend) - Uterm - N;

    Fields->AA.val[iPoint] = B1 * a01_times_invPow1_dt + B2 * a02_times_invPow2_dt - BM * divq * a01_times_invPow1_dt * (1.0 - BMblend) -
                             BM * (2.0 * DNA_POW2(divq) + detJ * MovingBoundary->Udot) * (1.0 - BMblend);

    /** RHS only involves the explicitly computed Laplacian **/
    Fields->RHS.val[iPoint] = -BL * Fields->dXI2_phi.val[iPoint] - BLr * Fields->dXI1_phi.val[iPoint];
  }

  return 0;
}

int TransformEqn_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties)
{
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
  DNA_FLOAT K = 2.0 * (FluidProperties->beta - 1.0 * NLflag) / pow2_c0;

  DNA_FLOAT detJ = Fields->Grid.detJacobi;

  for (int iPoint = 0; iPoint < RunOptions->NumericsFD.NPoints; iPoint++)
  {
    q = Fields->Grid.q[iPoint];

    geometryFactor = RunOptions->GeometricalDecay(Fields->Grid.x[iPoint]);

    u = Fields->BackgroundVelocity.val[iPoint];

    transformed_dXI1_phi = detJ * Fields->Old_dXI1_phi.o[0].val[iPoint];
    transformed_dt1 = Fields->dt1_phi.val[iPoint] + q * Fields->Old_dXI1_phi.o[0].val[iPoint];

    UL = 2.0 * u * detJ * q;

    A2 = 1.0 + K * u * transformed_dXI1_phi;

    B2 = A2 + K * transformed_dt1;

    AL = DNA_POW2(u) - pow2_c0 + 2.0 * u * transformed_dXI1_phi * NLflag + K * DNA_POW3(u) * transformed_dXI1_phi;

    BL = B2 * DNA_POW2(q) + AL * DNA_POW2(detJ) + UL;

    BLr = -pow2_c0 * geometryFactor * detJ;

    /** Assembling new RHS for the corrector step (involving the upated Laplacian) **/
    Fields->RHS.val[iPoint] = -BL * Fields->dXI2_phi.val[iPoint] - BLr * Fields->dXI1_phi.val[iPoint];
  }

  return 0;
}
