#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following funcitons prescribe the instantaneous position (R), velocity (U), and
acceleration (Udot) of the moving domain boundary.
---------------------------------------------------------**/

/** Constant velocity of the moving boundary **/
int BoundaryMotionLinear(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                         DNA_FLOAT time)
{
  MovingBoundary->R =
      MovingBoundary->R0 + RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude * (time - RunOptions->WaveExcitation.BoundaryMotionStartTime);
  MovingBoundary->U = RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude;
  MovingBoundary->Udot = 0.0;

  return 0;
}

/** Sinusoidal motion of the moving boundary **/
int BoundaryMotionOscillating(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                              DNA_FLOAT time)
{
  DNA_FLOAT Omega = (DNA_FLOAT) 2.0 * DNA_PI * RunOptions->WaveExcitation.MovingBoundaryFrequency;

  MovingBoundary->R = MovingBoundary->R0 + RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude / (Omega + DNA_SMALL) *
                                               DNA_SIN((time - RunOptions->WaveExcitation.BoundaryMotionStartTime) * Omega);
  MovingBoundary->U = RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude * DNA_COS((time - RunOptions->WaveExcitation.BoundaryMotionStartTime) * Omega);
  MovingBoundary->Udot = -RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude * (Omega + DNA_SMALL) *
                         DNA_SIN((time - RunOptions->WaveExcitation.BoundaryMotionStartTime) * Omega);

  return 0;
}

/**--------------------------------------------------------- 
Boundary motion that establishes a stationary background flow field with sonic
transition at a pre-defined sonic horizon radius rh (RunOptions->radius_horizon).
The flow field is based on the equation u0(r,t) = dotR*[R(t)/r]^k, where k=2
represents an incompressible flow field. The physically exact representation of
a stationary spherically symmetric flow field with smooth sonic transition
(promoted by simplified models of Bondi-Parker type flows) would require a
certain spatial distribution k(r). However, such a representation of the
compressible background flow field would also require to assume non-constant rho0
and c0 in the derivation of the governing wave equation. Therefore, k=const in the
present framework, so that the effect of a compressible flow field can at least
be mimicked in some close vicinity of the sonic horizon.
---------------------------------------------------------**/
int BoundaryMotionStationaryBlackHole(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                                      DNA_FLOAT time)
{
  RunOptions->statHorizon = 1;
  DNA_FLOAT k = RunOptions->U_backgroundScaling;
  DNA_FLOAT kplusone = k + 1.0;
  DNA_FLOAT one_by_kplusone = 1.0 / kplusone;
  DNA_FLOAT min_k_by_kplusone = -k / kplusone;
  DNA_FLOAT min_2kplusone_by_kplusone = -(2.0 * k + 1.0) / kplusone;
  DNA_FLOAT arg = DNA_POW(MovingBoundary->R0, kplusone) - kplusone * FluidProperties->c0 * DNA_POW(RunOptions->radius_horizon, k) * time;

  MovingBoundary->R = DNA_POW(arg, one_by_kplusone);
  MovingBoundary->U = -FluidProperties->c0 * DNA_POW(RunOptions->radius_horizon, k) * DNA_POW(arg, min_k_by_kplusone);
  MovingBoundary->Udot = -k * DNA_POW2(FluidProperties->c0) * DNA_POW(RunOptions->radius_horizon, 2.0 * k) * DNA_POW(arg, min_2kplusone_by_kplusone);

  return 0;
}

int BoundaryMotionStationaryWhiteHole(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                                      DNA_FLOAT time)
{
  RunOptions->statHorizon = 1;
  DNA_FLOAT k = RunOptions->U_backgroundScaling;
  DNA_FLOAT kplusone = k + 1.0;
  DNA_FLOAT one_by_kplusone = 1.0 / kplusone;
  DNA_FLOAT min_k_by_kplusone = -k / kplusone;
  DNA_FLOAT min_2kplusone_by_kplusone = -(2.0 * k + 1.0) / kplusone;
  DNA_FLOAT arg = DNA_POW(MovingBoundary->R0, kplusone) + kplusone * FluidProperties->c0 * DNA_POW(RunOptions->radius_horizon, k) * time;

  MovingBoundary->R = DNA_POW(arg, one_by_kplusone);
  MovingBoundary->U = FluidProperties->c0 * DNA_POW(RunOptions->radius_horizon, k) * DNA_POW(arg, min_k_by_kplusone);
  MovingBoundary->Udot = -k * DNA_POW2(FluidProperties->c0) * DNA_POW(RunOptions->radius_horizon, 2.0 * k) * DNA_POW(arg, min_2kplusone_by_kplusone);

  return 0;
}
