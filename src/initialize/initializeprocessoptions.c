#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
In this functions, all function pointers as declared in <DNA.h> are set
based on the settings specified in the options file (see ioreadoptionsfile.c).
---------------------------------------------------------**/

int InitializeProcessOptions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{ 
  if (RunOptions->NumericsFD.BCEast == 0)
  {
    RunOptions->FDBC_East = &FDBC_Mur_dummy;
  }
  if (RunOptions->NumericsFD.BCEast == 1)
  {
    if (MovingBoundary->R0 < Fields->Grid.xfix) RunOptions->FDBC_East = &FDBC_Mur_East;
    else RunOptions->FDBC_East = &FDBC_Mur_East_inv;
  }
  if (RunOptions->NumericsFD.BCWest == 0)
  {
    RunOptions->FDBC_West = &FDBC_Mur_dummy;
  }
  if (RunOptions->NumericsFD.BCWest == 1)
  {
    if (MovingBoundary->R0 < Fields->Grid.xfix) RunOptions->FDBC_West = &FDBC_Mur_West;
    else RunOptions->FDBC_West = &FDBC_Mur_West_inv;
  }
  
  if (RunOptions->WaveExcitation.GaussEnvelopeOption == 0)
  {
    RunOptions->GaussEnvelope = &WaveExcitationGaussEnvelopeNone;
  }
  if (RunOptions->WaveExcitation.GaussEnvelopeOption == 1)
  {
    RunOptions->GaussEnvelope = &WaveExcitationGaussEnvelope;
  }

  if (RunOptions->WaveExcitation.ExcitationFunction == 0)
  {
    RunOptions->PeriodicPressureExcitation = &WaveExcitationSineModulated;
  }
  if (RunOptions->WaveExcitation.ExcitationFunction == 1)
  {
    RunOptions->PeriodicPressureExcitation = &WaveExcitationConst;
    /** Make sure that this is not be overwritten by the setting for the GaussEnvelopeOption **/
    RunOptions->GaussEnvelope = &WaveExcitationGaussEnvelope;
  }

  if (RunOptions->BoundaryMotionType == 0)
  {
    RunOptions->BoundaryMotion = &BoundaryMotionLinear;
  }
  if (RunOptions->BoundaryMotionType == 1)
  {
    RunOptions->BoundaryMotion = &BoundaryMotionOscillating;
  }
  if (RunOptions->BoundaryMotionType == 2)
  {
    RunOptions->BoundaryMotion = &BoundaryMotionStationaryBlackHole;
  }
  if (RunOptions->BoundaryMotionType == 3)
  {
    RunOptions->BoundaryMotion = &BoundaryMotionStationaryWhiteHole;
  }
  if (RunOptions->BoundaryMotionType == 2 || RunOptions->BoundaryMotionType == 3)
  {
    if (RunOptions->GravitationalPotentialOption == 0)
    {
      RunOptions->WaveCalcGravitationalPotential = &BackgroundFlowGravitationalPotential_dummy;
    }
    if (RunOptions->GravitationalPotentialOption == 1)
    {
      RunOptions->WaveCalcGravitationalPotential = &BackgroundFlowGravitationalPotential;
    }
  }
  
  if (RunOptions->backgroundMotionMode == 0)
  {
    RunOptions->WaveBackgroundVelocityAtWall = &BackgroundFlowVelocityAtWall_decoupled;
    RunOptions->WaveBackgroundAccelerationAtWall = &BackgroundFlowAccelerationAtWall_decoupled;
  }
  if (RunOptions->backgroundMotionMode == 1)
  {
    RunOptions->WaveBackgroundVelocityAtWall = &BackgroundFlowVelocityAtWall_coupledToWall;
    RunOptions->WaveBackgroundAccelerationAtWall = &BackgroundFlowAccelerationAtWall_coupledToWall;
  }
  
  if (RunOptions->dimension == 1)
  {
    if (RunOptions->backgroundMotionMode == 0){RunOptions->WaveBackgroundMotion = &BackgroundFlowMotion_const;}
    if (RunOptions->backgroundMotionMode == 1){RunOptions->WaveBackgroundMotion = &BackgroundFlowMotion_Cartesian;}
    RunOptions->GeometricalDecay = &TransformNoGeometricalDecay;
  }
  if (RunOptions->dimension == 3)
  {
    RunOptions->WaveBackgroundMotion = &BackgroundFlowMotion_spherical;
    RunOptions->GeometricalDecay = &TransformSphericalDecay;
  }

  if(RunOptions->Probes.nSamplePoints > 0)
  {
    RunOptions->IOWriteProbesOption = &IOWriteProbes;
    RunOptions->IOUpdateProbesOption = &IOUpdateProbes;
  }
  else
  {
    RunOptions->IOWriteProbesOption = &IOWriteProbes_dummy;
    RunOptions->IOUpdateProbesOption = &IOUpdateProbes_dummy;
  }
  if(RunOptions->BoundaryMotionType == 2 || RunOptions->BoundaryMotionType == 3)
  {
    RunOptions->IOWriteStatHorizonOption = &IOWriteStatHorizon;
  }
  else
  {
    RunOptions->IOWriteStatHorizonOption = &IOWriteStatHorizon_dummy;
  }

  return 0;
}