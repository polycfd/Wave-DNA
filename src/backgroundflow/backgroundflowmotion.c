#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions compute the background flow field u0 and the required 
derivatives thereof based on the geometry of the flow problem and the 
instantaneous position, velocity, and acceleration of the moving domain boundary.
---------------------------------------------------------**/

/** Spatially and temporarily constant background flow field **/
int BackgroundFlowMotion_const(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.BackgroundVelocity[iPoint] = MovingBoundary->U_backgroundAtWall;
    Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] = 0.0;
    Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] = 0.0;
    Fields->BackgroundFlowField.dt1material_BackgroundVelocity[iPoint] = 0.0;
  }

  return 0;
}

/** Transient and spherically symmetric background flow field **/
int BackgroundFlowMotion_spherical(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  DNA_FLOAT k = RunOptions->U_backgroundScaling;

  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.BackgroundVelocity[iPoint] = DNA_POW(MovingBoundary->R / (Fields->Grid.x[iPoint]), k) * MovingBoundary->U_backgroundAtWall;
  }

  for (register int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] = -k * Fields->BackgroundFlowField.BackgroundVelocity[iPoint] / (Fields->Grid.x[iPoint]);

    Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] =
        DNA_POW(MovingBoundary->R / (Fields->Grid.x[iPoint]), k) * MovingBoundary->Udot +
        k * MovingBoundary->U_backgroundAtWall / MovingBoundary->R * Fields->BackgroundFlowField.BackgroundVelocity[iPoint];
  }

  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.dt1material_BackgroundVelocity[iPoint] =
        Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] +
        Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] * Fields->BackgroundFlowField.BackgroundVelocity[iPoint];
  }

  return 0;
}

/** Transient but spatially uniform background flow field (1d Cartesian coordinates) **/
int BackgroundFlowMotion_Cartesian(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.BackgroundVelocity[iPoint] = MovingBoundary->U_backgroundAtWall;
  }

  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] = 0.0;

    Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] = MovingBoundary->Udot;
  }

  for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BackgroundFlowField.dt1material_BackgroundVelocity[iPoint] =
        Fields->BackgroundFlowField.dt1_BackgroundVelocity[iPoint] +
        Fields->BackgroundFlowField.GradBackgroundVelocity[iPoint] * Fields->BackgroundFlowField.BackgroundVelocity[iPoint];
  }

  return 0;
}

/**--------------------------------------------------------- 
The background flow velocity and acceleration at the moving domain boundary are set 
equal to the boundary velocity and acceleration, respectively, if the background 
flow field is chosen to be coupled to the boundary motion. In that case, the moving 
boundary mimics a solid boundary that displaces the surrounding fluid. If the flow
is chosen to be decoupled, the domain boundary is a virtual boundary moving relative
to the fluid. 
---------------------------------------------------------**/

DNA_FLOAT BackgroundFlowVelocityAtWall_coupledToWall(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  return MovingBoundary->U;
}

DNA_FLOAT BackgroundFlowVelocityAtWall_decoupled(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time) { return 0.0; }

DNA_FLOAT BackgroundFlowAccelerationAtWall_coupledToWall(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  return MovingBoundary->Udot;
}

DNA_FLOAT BackgroundFlowAccelerationAtWall_decoupled(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time)
{
  return 0.0;
}
