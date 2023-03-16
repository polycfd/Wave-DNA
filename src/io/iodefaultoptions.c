#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Setting the default options to enable a default run even for an empty
options file.
---------------------------------------------------------**/

int IODefaultOptions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                     struct DNA_FluidProperties *FluidProperties)
{
  /** Fluid properties **/
  FluidProperties->c0 = 1500.0;
  FluidProperties->rho0 = 1000.0;
  RunOptions->LagrangianDensityFlag = 0.0;
  FluidProperties->beta = 0.0;

  /** Wave excitation **/
  RunOptions->iExcitation = 0;
  RunOptions->WaveExcitation.PressureAmplitude = 1e3;
  RunOptions->WaveExcitation.ExcitationFrequency = 1e5;
  RunOptions->WaveExcitation.PhaseShiftAngle = 0;
  RunOptions->WaveExcitation.PowerCoeff = 1;
  RunOptions->WaveExcitation.EnvelopedPeriod = 1;
  RunOptions->WaveExcitation.GaussShapeCoeff = 1.0;

  /** Physical domain **/
  // The default length of the domain comprises 10 emitted wavelengths
  MovingBoundary->R0 = 0.0;
  Fields->Grid.xfix = 10.0 * FluidProperties->c0 / RunOptions->WaveExcitation.ExcitationFrequency;
  Fields->Grid.xmov = MovingBoundary->R0;

  /** Numerical resolution **/
  // 100 grid points per wavelength
  RunOptions->NumericsFD.NPoints =
      100 * (int) ((Fields->Grid.xfix - Fields->Grid.xmov) / (FluidProperties->c0 / RunOptions->WaveExcitation.ExcitationFrequency)) + 1;
  RunOptions->NumericsFD.dx = (Fields->Grid.xfix - Fields->Grid.xmov) / ((DNA_FLOAT) RunOptions->NumericsFD.NPoints - 1);
  // CFL0 = 0.1 per default
  RunOptions->NumericsFD.dt = RunOptions->NumericsFD.dx * 0.1 / (FluidProperties->c0);

  /** Physical simulation time **/
  RunOptions->tStart = 0.0;
  // Per default, 10 emitted wave periods are simulated
  RunOptions->tEnd = 10.0 / RunOptions->WaveExcitation.ExcitationFrequency;
  RunOptions->WaveExcitation.ExcitationStart = RunOptions->tStart;
  RunOptions->WaveExcitation.ExcitationEnd = RunOptions->tEnd;
  RunOptions->WaveExcitation.BoundaryMotionStartTime = RunOptions->tStart;
  RunOptions->WaveExcitation.BoundaryMotionEndTime = RunOptions->tEnd;

  /** Boundary motion **/
  MovingBoundary->R0 = 0;
  MovingBoundary->U_backgroundAtWall = 0;
  MovingBoundary->Udot_backgroundAtWall = 0;
  RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude = 0;
  RunOptions->WaveExcitation.MovingBoundaryFrequency = 0.1 * RunOptions->WaveExcitation.ExcitationFrequency;

  /** Background flow **/
  RunOptions->statHorizon = 0;
  RunOptions->radius_horizon = 0.9 * Fields->Grid.xmov;
  RunOptions->U_backgroundScaling = 2.0;
  RunOptions->GravitationalPotentialOption = 0;

  /** Finite difference **/
  RunOptions->NumericsFD.nCorrectors = 1;
  RunOptions->NumericsFD.correctorWeight = 0.1;

  /** Results **/
  // The default write frequency is set in such a way that the results are written only once
  // at the end of the simulation
  RunOptions->writeFrequency = (int) ((RunOptions->tEnd - RunOptions->tStart) / RunOptions->NumericsFD.dt);
  RunOptions->Probes.nSamplePoints = 0;

  /** The following options are needed to set the function pointers (process options) **/
  RunOptions->NumericsFD.BCEast = 1;
  RunOptions->NumericsFD.BCWest = 1;
  RunOptions->WaveExcitation.GaussEnvelopeOption = 0;
  RunOptions->WaveExcitation.ExcitationFunction = 0;
  RunOptions->backgroundMotionMode = 0;
  RunOptions->dimension = 1;
  RunOptions->BoundaryMotionType = 0;

  RunOptions->CalcGravitationalPotential = &BackgroundFlowGravitationalPotential_dummy;

  return 0;
}
