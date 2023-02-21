#include "DNA.h"
#include "DNA-functions.h"

/** Called in main.c if only wave is to be solved **/
int SolveTimeLoop(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties)
{
  RunOptions->t = RunOptions->tStart;
  RunOptions->tEnd = RunOptions->tEnd;
  RunOptions->NumericsFD.dtNumber = 0;

  while (RunOptions->t < (RunOptions->tEnd + DNA_EPS))
  {
    /** Calling the function to prescribe the boundary motion **/
    if(
           (RunOptions->t >= RunOptions->WaveExcitation.BoundaryMotionStartTime)
        && (RunOptions->t <= RunOptions->WaveExcitation.BoundaryMotionEndTime) 
      )
    {
      (*RunOptions->BoundaryMotion)(RunOptions, MovingBoundary, FluidProperties, RunOptions->t);
    }
    else
    {
      // MovingBoundary->R remains unchanged
      MovingBoundary->U = 0.0;
      MovingBoundary->Udot = 0.0;
    }

    /** Value by which the pressure excitation signal is Gauss convoluted **/
    RunOptions->WaveExcitation.GaussConvolutionFactor = (*RunOptions->GaussEnvelope)(
        RunOptions->WaveExcitation.EnvelopedPeriod, RunOptions->WaveExcitation.ExcitationFrequency, RunOptions->t, RunOptions->WaveExcitation.GaussShapeCoeff);

    /** setting the (periodic) excitation pressure as given by a predefined function **/
    RunOptions->pExcitation = (*RunOptions->PeriodicPressureExcitation)(RunOptions);

    Solve(RunOptions, Fields, MovingBoundary, FluidProperties);
    
    RunOptions->t += RunOptions->NumericsFD.dt;
    ++(RunOptions->NumericsFD.dtNumber);
  }

  return 0;
}




/** This function is always called when wave is to be solved **/
int Solve(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties)
{
  int corrNo;

  /** Sum over previous function values times corresponding FD coefficients **/
  SolveSumAy_dt(RunOptions, Fields);

   /** Updating the position of the excitation boundary and the uniform grid (metrics) **/
  Fields->Grid.xmov = MovingBoundary->R;
  GridPhysicalDomain(&(RunOptions->NumericsFD), &(Fields->Grid));
  GridMotion(RunOptions, Fields, MovingBoundary);

  /** Background motion at the wall and in the entire field **/
  MovingBoundary->U_backgroundAtWall = (*RunOptions->WaveBackgroundVelocityAtWall)(RunOptions, MovingBoundary, RunOptions->t);
  MovingBoundary->Udot_backgroundAtWall = (*RunOptions->WaveBackgroundAccelerationAtWall)(RunOptions, MovingBoundary, RunOptions->t);
  (*RunOptions->WaveBackgroundMotion)(RunOptions, Fields, MovingBoundary, RunOptions->t);

  //ToDo: check appropriate position of this
  SolveCalcqphi(RunOptions, Fields);
  
  /** Predictior step: coordinate transformations, explicit solution, storing results of initial guess **/
  TransformEqn_Predictor(RunOptions, Fields, MovingBoundary, FluidProperties);
  SolveExplicit_Predictor(RunOptions, Fields);
  SolveStoreInitialGuess(RunOptions, Fields);

  /** Loop over the corrector steps **/
  for (corrNo = 0; corrNo < RunOptions->NumericsFD.nCorrectors; corrNo++)
  {
    SolveCorrectLaplacian(RunOptions, Fields);
    TransformEqn_Corrector(RunOptions, Fields, FluidProperties);
    SolveExplicit_Corrector(RunOptions, Fields);
  }  

  /** Updating explcit spatial/time derivatives **/
  SolveExplicitDerivatives_dx(RunOptions, Fields);
  SolveExplicitDerivatives_dt(RunOptions, Fields);

  BoundaryConditions(RunOptions, Fields, FluidProperties);

  TransformPotentialFieldToPressureField(RunOptions, Fields, FluidProperties);
  WaveExcitationSetPressure(RunOptions, Fields, RunOptions->t);
  TransformExcitationPressureToExcitationPotential(RunOptions, Fields, FluidProperties);
  //RunOptions->phiField.val[RunOptions->iExcitation] = RunOptions->pExcitation;

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
