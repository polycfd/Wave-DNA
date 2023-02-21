#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>

#include "DNA-constants.h"

#ifndef DNA_H_
#define DNA_H_

/**--------------------------------------------------------- 
Declaration of all structures and (function) pointers. 
<DNA.h> must always be included.
---------------------------------------------------------**/

/** Pre-declaration of all structs
(needed as there are some nested structs) **/
struct DNA_FDCoeffs;
struct DNA_Grid;
struct DNA_Probes;
struct DNA_ScalarField;
struct DNA_OldScalarFields;
struct DNA_MovingBoundary;
struct DNA_RunOptions;
struct DNA_Fields;
struct DNA_FluidProperties;
struct DNA_NumericsFD;

/** Finite difference coefficients **/
struct DNA_FDCoeffs
{
  struct adt_order1
  {
    DNA_FLOAT dt1[2];
    DNA_FLOAT dt2[3];
  } adt_order1;

  struct adx_order2
  {
    DNA_FLOAT dx1[3];
    DNA_FLOAT dx2[3];
  } adx_order2;
};

/** Physical and computational grids and metrics of the coordinate transformation **/
struct DNA_Grid
{
  DNA_FLOAT *x;
  DNA_FLOAT *XI;
  DNA_FLOAT wallDisplacement;
  DNA_FLOAT xmov;
  DNA_FLOAT xfix;
  DNA_FLOAT *q;
  DNA_FLOAT *dtq;
  DNA_FLOAT divq;
  DNA_FLOAT detJacobi;
  DNA_FLOAT dxidetJacobi;
};

/** Local time-signals of the acoustic pressure can be taken at sample points **/
struct DNA_Probes
{
  int nSamplePoints;
  int *SampleIDs;
  int *TimeIDs;
  DNA_FLOAT *SamplePoints;
  DNA_FLOAT *SampleTime;
  DNA_FLOAT **SamplePressure;
};

/** All fields (pressure or velocity field etc) are instances of this struct **/
struct DNA_ScalarField
{
  DNA_FLOAT *val;
};

/** Similar to the above, but contains information of the old time steps **/
struct DNA_OldScalarFields
{
  struct DNA_ScalarField o[2];
};

/** Parameters of the acoustic pressure wave excitation **/
struct DNA_WaveExcitation
{
  int ExcitationFunction;
  int GaussEnvelopeOption;
  DNA_FLOAT ExcitationStart;
  DNA_FLOAT ExcitationEnd;
  DNA_FLOAT PressureAmplitude;
  DNA_FLOAT ExcitationFrequency;
  DNA_FLOAT PhaseShiftAngle;
  DNA_FLOAT PowerCoeff;
  DNA_FLOAT EnvelopedPeriod;
  DNA_FLOAT GaussShapeCoeff;
  DNA_FLOAT GaussConvolutionFactor;
  DNA_FLOAT MovingBoundaryVelocityAmplitude;
  DNA_FLOAT MovingBoundaryFrequency;
  DNA_FLOAT BoundaryMotionStartTime;
  DNA_FLOAT BoundaryMotionEndTime;
};

/** Parameters specifiying the numerical resolution, the boundary conditions,
and the settings of the predictor-corrector method **/
struct DNA_NumericsFD
{
  int dtNumber;
  DNA_FLOAT dt;
  DNA_FLOAT dx;
  DNA_FLOAT dXI;
  int NPoints;
  int BCWest;
  int BCEast;
  int nCorrectors;
  DNA_FLOAT correctorWeight; 

  struct DNA_FDCoeffs FDCoeffs;
};

/** Run options of the simulation **/
struct DNA_RunOptions
{ 
  int OnScreenIOPriority;  // Determines how much stuff is printed in the terminal (low = pretty much all, high = only essentials)
  char OptionsDir[DNA_STRINGLENGTH];  // Directory path of the options file

  int dimension;
  int backgroundMotionMode;
  int BoundaryMotionType;
  int statHorizon;
  int GravitationalPotentialOption;
  
  DNA_FLOAT tStart;
  DNA_FLOAT tEnd;
  DNA_FLOAT t;
  
  int LagrangianDensityFlag;
  DNA_FLOAT radius_horizon;
  DNA_FLOAT U_backgroundScaling;

  int writeCount;
  int writeFrequency;
  FILE *WaveResults;
  int (*IOWriteProbesOption)(struct DNA_RunOptions *RunOptions);
  int (*IOUpdateProbesOption)(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
  int (*IOWriteStatHorizonOption)(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);

  struct DNA_Probes Probes;
  struct DNA_NumericsFD NumericsFD;

  /** Wave excitation **/
  int iExcitation;
  DNA_FLOAT pExcitation;
  struct DNA_WaveExcitation WaveExcitation;

  // Pointers to the pressure excitation function and the Gauss envelope function
  DNA_FLOAT (*PeriodicPressureExcitation)(struct DNA_RunOptions *RunOptions);
  DNA_FLOAT (*GaussEnvelope)(DNA_FLOAT periodNo, DNA_FLOAT frequ, DNA_FLOAT time, DNA_FLOAT shapeCoeff);

  /** Pointers to the boundary condition functions **/
  int (*FDBC_East)(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, struct DNA_ScalarField *NewField, struct DNA_ScalarField *OldField);
  int (*FDBC_West)(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, struct DNA_ScalarField *NewField, struct DNA_ScalarField *OldField);

  // Pointer to the function describing the decay/amplification due to the geometry of the problem
  double (*GeometricalDecay)(DNA_FLOAT x);

  // Pointers to the boundary and background flow motion functions
  int (*BoundaryMotion)(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT time);
  int (*WaveBackgroundMotion)(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
  double (*WaveBackgroundVelocityAtWall)(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
  double (*WaveBackgroundAccelerationAtWall)(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
  double (*WaveCalcGravitationalPotential)(struct DNA_RunOptions *RunOptions, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT x);
};

/** Fluid properties **/
struct DNA_FluidProperties
{
  DNA_FLOAT c0;
  DNA_FLOAT rho0;
  DNA_FLOAT beta;
};

/** Parameters describing the boundary motion
(left/West boundary of the physical domain) **/
struct DNA_MovingBoundary
{
  DNA_FLOAT R0;
  DNA_FLOAT R;
  DNA_FLOAT U;
  DNA_FLOAT Udot;
  DNA_FLOAT U_backgroundAtWall;
  DNA_FLOAT Udot_backgroundAtWall;
};

/** Struct containing all fields **/
struct DNA_Fields
{
  int sizeof_OldScalarFields;

  struct DNA_ScalarField phi;
  struct DNA_ScalarField qphi;
  struct DNA_ScalarField dXI1_phi;
  struct DNA_ScalarField dXI2_phi;
  struct DNA_ScalarField dXI1_qphi;
  struct DNA_ScalarField dt1_phi;
  struct DNA_ScalarField dt2_phi;
  struct DNA_ScalarField dt1_dXI1_phi;
  struct DNA_ScalarField dt1_dXI1_qphi;
  struct DNA_ScalarField PressureField;
  
  struct DNA_ScalarField sum_aphi_dt1;
  struct DNA_ScalarField sum_aphi_dt2;
  struct DNA_ScalarField sum_adXI1_phi_dt1;
  struct DNA_ScalarField sum_adXI1_qphi_dt1;
  
  struct DNA_ScalarField phi1_initGuess;
  struct DNA_ScalarField RHS;
  struct DNA_ScalarField AA;
  struct DNA_ScalarField BB;
  struct DNA_ScalarField BBbyAA;

  struct DNA_ScalarField BackgroundVelocity;
  struct DNA_ScalarField GradBackgroundVelocity;
  struct DNA_ScalarField dt1_BackgroundVelocity;
  struct DNA_ScalarField dt1material_BackgroundVelocity;

  struct DNA_OldScalarFields Old_phi;
  struct DNA_OldScalarFields Old_dXI1_phi;
  struct DNA_OldScalarFields Old_dXI2_phi;
  struct DNA_OldScalarFields Old_dXI1_qphi;

  struct DNA_Grid Grid;
};

#endif /* DNA_H_ */
