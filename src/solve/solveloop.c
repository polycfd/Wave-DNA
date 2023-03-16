#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The fucntion <SolveTimeLoop> is called in main.c and represents the time loop of the
simulation. In <SolveTimeLoop>, the function <Solve> is called, in which the numerical
algorithm is implemented. 
---------------------------------------------------------**/

int SolveTimeLoop(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                  struct DNA_FluidProperties *FluidProperties)
{
  RunOptions->t = RunOptions->tStart;
  RunOptions->tEnd = RunOptions->tEnd;
  RunOptions->NumericsFD.dtNumber = 0;

  int progress = 0;
  IOProgressInitial();

  while (RunOptions->t < (RunOptions->tEnd + DNA_EPS))
  {
    /** Calling the function to prescribe the boundary motion **/
    if ((RunOptions->t > (RunOptions->WaveExcitation.BoundaryMotionStartTime - DNA_SMALL)) &&
        (RunOptions->t < (RunOptions->WaveExcitation.BoundaryMotionEndTime + DNA_SMALL)))
    {
      RunOptions->BoundaryMotion(RunOptions, MovingBoundary, FluidProperties, RunOptions->t);
    }
    else
    {
      // MovingBoundary->R remains unchanged
      MovingBoundary->U = 0.0;
      MovingBoundary->Udot = 0.0;
    }

    /** Value by which the pressure excitation signal is Gauss convoluted **/
    RunOptions->WaveExcitation.GaussConvolutionFactor = RunOptions->GaussEnvelope(
        RunOptions->WaveExcitation.EnvelopedPeriod, RunOptions->WaveExcitation.ExcitationFrequency, RunOptions->t, RunOptions->WaveExcitation.GaussShapeCoeff);

    /** setting the (periodic) excitation pressure as given by a predefined function **/
    RunOptions->pExcitation = RunOptions->PeriodicPressureExcitation(RunOptions);

    Solve(RunOptions, Fields, MovingBoundary, FluidProperties);

    RunOptions->t += RunOptions->NumericsFD.dt;
    ++(RunOptions->NumericsFD.dtNumber);

    IOProgressUpdate(&progress, RunOptions->t, RunOptions->tEnd);
  }

  IOProgressFinal();

  return 0;
}

int Solve(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties)
{
  /** Sum over previous function values times corresponding FD coefficients **/
  SolveSumAy_dt(RunOptions, Fields);

  /** Updating the position of the excitation boundary and the uniform grid (metrics) **/
  Fields->Grid.xmov = MovingBoundary->R;
  GridPhysicalDomain(&RunOptions->NumericsFD, &Fields->Grid);
  GridMotion(RunOptions, Fields, MovingBoundary);

  /** Background motion at the wall and in the entire field **/
  MovingBoundary->U_backgroundAtWall = RunOptions->BackgroundVelocityAtWall(RunOptions, MovingBoundary, RunOptions->t);
  MovingBoundary->Udot_backgroundAtWall = RunOptions->BackgroundAccelerationAtWall(RunOptions, MovingBoundary, RunOptions->t);
  RunOptions->BackgroundMotion(RunOptions, Fields, MovingBoundary, RunOptions->t);

  /** The term qphi represents the product of the acoustic perturbation potential phi
  and the velocity q = dXI/dt of the grid points in the computational domain as viewed
  from a moving point in the physical domain (also taking the Jacobain into ccoutn) **/
  SolveCalcqphi(RunOptions, Fields);

  /** Predictior step: coordinate transformations, explicit solution, storing results of initial guess **/
  TransformEqn_Predictor(RunOptions, Fields, MovingBoundary, FluidProperties);
  SolveExplicit_Predictor(RunOptions, Fields);
  SolveStoreInitialGuess(RunOptions, Fields);

  /** Loop over the corrector steps **/
  for (int corrNo = 0; corrNo < RunOptions->NumericsFD.nCorrectors; corrNo++)
  {
    SolveCorrectLaplacian(RunOptions, Fields);
    TransformEqn_Corrector(RunOptions, Fields, FluidProperties);
    SolveExplicit_Corrector(RunOptions, Fields);
  }

  /** Updating explcit spatial/time derivatives **/
  SolveExplicitDerivatives_dx(RunOptions, Fields);
  SolveExplicitDerivatives_dt(RunOptions, Fields);

  BoundaryConditions(RunOptions, Fields, FluidProperties);

  /** The following block converts the acoustic perturbation potential field
  (primary solution of the numerical procedure), sets the new excitation pressure
  at the emission node and subsequently converts the excitation pressure into the
  excitation perturbation potential (typically a boundary condition). **/
  TransformPotentialFieldToPressureField(RunOptions, Fields, FluidProperties);
  TransformExcitationPressureToExcitationPotential(RunOptions, Fields, FluidProperties);

  /** Updating the old fields required for time integration **/
  SolveUpdateOldFields(RunOptions, Fields);

  /** Write results to disc and update write counter **/
  if (RunOptions->writeCount == 0)
  {
    IOWriteToDisc(RunOptions, Fields, MovingBoundary);
    if (RunOptions->NumericsFD.dtNumber > 0)
    {
      RunOptions->IOWriteProbesOption(RunOptions);
    }
  }

  RunOptions->IOUpdateProbesOption(RunOptions->writeCount, RunOptions, Fields);
  RunOptions->IOWriteStatHorizonOption(RunOptions, Fields, MovingBoundary);
  RunOptions->writeCount++;

  if (RunOptions->writeCount == RunOptions->writeFrequency)
  {
    RunOptions->writeCount = 0;
  }

  return 0;
}
