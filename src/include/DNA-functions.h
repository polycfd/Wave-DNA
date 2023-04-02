#include "DNA.h"

#ifndef DNAFUNCTIONS_H_
#define DNAFUNCTIONS_H_

/**--------------------------------------------------------- 
Prototypes of all functions that are used in Wave-DNA. The header 
<DNA-functions.h> must be included to enable function calls. Some
of the functions are called by reference. The corresponding function
pointers are declared in <DNA.h> and set in initializeprocessoptions.c 
based on the settings specified in the options file.
---------------------------------------------------------**/

/** Initializion of the simulation **/
int InitializeSimulation(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);
int InitializeProcessOptions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);
int InitializeScalarField(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *ScalarField, DNA_FLOAT val);
int InitializeOldPhiFields(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT **OldPhiFields, DNA_FLOAT val);

/** Wave excitation **/
DNA_FLOAT WaveExcitationSineModulated(struct DNA_RunOptions *RunOptions);
DNA_FLOAT WaveExcitationConst(struct DNA_RunOptions *RunOptions);
DNA_FLOAT WaveExcitationGaussEnvelope(DNA_FLOAT periodNo, DNA_FLOAT frequ, DNA_FLOAT time, DNA_FLOAT shapeCoeff);
DNA_FLOAT WaveExcitationGaussEnvelopeNone(DNA_FLOAT periodNo, DNA_FLOAT frequ, DNA_FLOAT time, DNA_FLOAT shapeCoeff);

/** Transform between acoustic potential and pressure **/
int TransformExcitationPressureToExcitationPotential(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties);
int TransformPotentialFieldToPressureField(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties);

/** Transformation between the physical and the computational domain and geometrical decay in spherical symmetry **/
int TransformEqn_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                           struct DNA_FluidProperties *FluidProperties);
int TransformEqn_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties);
DNA_FLOAT TransformSphericalDecay(DNA_FLOAT x);
DNA_FLOAT TransformNoGeometricalDecay(DNA_FLOAT x);

/** Functions to compute the properties of the background flow field **/
int BackgroundFlowMotion_const(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
int BackgroundFlowMotion_spherical(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
int BackgroundFlowMotion_Cartesian(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
DNA_FLOAT BackgroundFlowVelocityAtMovingBoundary_coupledToWall(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
DNA_FLOAT BackgroundFlowAccelerationAtMovingBoundary_decoupled(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
DNA_FLOAT BackgroundFlowAccelerationAtMovingBoundary_coupledToWall(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary,
                                                                   DNA_FLOAT time);
DNA_FLOAT BackgroundFlowVelocityAtMovingBoundary_decoupled(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, DNA_FLOAT time);
DNA_FLOAT BackgroundFlowGravitationalPotential(struct DNA_RunOptions *RunOptions, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT x);
DNA_FLOAT BackgroundFlowGravitationalPotential_dummy(struct DNA_RunOptions *RunOptions, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT x);

/** Functions related to the numerical solution algorithm **/
int Solve(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties);
int SolveTimeLoop(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                  struct DNA_FluidProperties *FluidProperties);
int SolveUpdateOldFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveCalcqphi(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveSumFDcoeffTimesPhi_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveExplicitDerivatives_dx(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveExplicitDerivatives_dt(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveExplicit_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveExplicit_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveCorrectLaplacian(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int SolveStoreInitialGuess(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);

/** Functions related to the finite difference discretization **/
int FDFiniteDifferenceCoeffs(struct DNA_NumericsFD *NumericsFD);
int FDGradient(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *Y, DNA_FLOAT *gradField);
int FDLaplacian(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *Y, DNA_FLOAT *gradField);
int FDdt1(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT *Y, DNA_FLOAT *dt1Field);
int FDdt2(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT *Y, DNA_FLOAT *dt2Field);
int FDSumFDcoeffTimesPhi_dt1(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT **OldPhiFields);
int FDSumFDcoeffTimesPhi_dt2(struct DNA_NumericsFD *NumericsFD, DNA_FLOAT *SumFDcoeffTimesPhi, DNA_FLOAT **OldPhiFields);
int FDBC_Mur_East(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT *NewField,
                  DNA_FLOAT *OldField);
int FDBC_Mur_West(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT *NewField,
                  DNA_FLOAT *OldField);
int FDBC_Mur_West_inv(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT *NewField,
                      DNA_FLOAT *OldField);
int FDBC_Mur_East_inv(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT *NewField,
                      DNA_FLOAT *OldField);
int FDBC_Mur_dummy(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties, DNA_FLOAT *NewField,
                   DNA_FLOAT *OldField);

/** Functions to set the boundary conditions and to describe the motion of the moving domain boundary **/
int BoundaryConditions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_FluidProperties *FluidProperties);
int BoundaryMotionLinear(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                         DNA_FLOAT time);
int BoundaryMotionOscillating(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                              DNA_FLOAT time);
int BoundaryMotionStationaryBlackHole(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                                      DNA_FLOAT time);
int BoundaryMotionStationaryWhiteHole(struct DNA_RunOptions *RunOptions, struct DNA_MovingBoundary *MovingBoundary, struct DNA_FluidProperties *FluidProperties,
                                      DNA_FLOAT time);

/** Functions to create the computational and physical grids and to compute the metrics of the coordinate transformation **/
int GridPhysicalDomain(struct DNA_NumericsFD *NumericsFD, struct DNA_Grid *Grid);
int GridComputationalDomain(struct DNA_NumericsFD *NumericsFD, struct DNA_Grid *Grid);
int GridMotion(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);

/** Field Memory Allocation **/
int MemoryAllocSamplePoints(int NPoints, struct DNA_Probes *Probes);
int MemoryAllocSampleVectors(int writeFrequency, struct DNA_Probes *Probes);
int MemoryAllocGrid(int NPoints, struct DNA_Grid *Grid);
int MemoryAllocBackgroundFlowField(int NPoints, struct DNA_BackgroundFlowField *BackgroundFlowField);
int MemoryAllocPhiField(int NPoints, struct DNA_PhiField *PhiField);
int MemoryAllocOldPhiField(int NPoints, struct DNA_OldPhiField *OldPhiField);
int MemoryAllocFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int MemoryFreeSampleVectors(int writeFrequency, struct DNA_Probes *Probes);
int MemoryFreeFields(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);

/** In put and output on screen and writing to disk **/
int IODefaultOptions(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                     struct DNA_FluidProperties *FluidProperties);
int IOErrorOnScreen(int num, char *message);
int IOWelcomeScreen();
int IOExitScreen(double totaltime);
int IOReadOneOption(FILE *OptionsFile, char *stuk);
int IOReadOptionsFile(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                      struct DNA_FluidProperties *FluidProperties);
int IOReadCommandLineOptions(int argc, char **args, struct DNA_RunOptions *RunOptions);
int IOLineGet(char *ssring, FILE *fp);
int IOLineGetSemi(char *ssring, FILE *fp);
int IOLineGetSkip(char *ssring, FILE *fp);
int IOWriteHeader(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int IOWriteToDisc(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);
int IOWriteStatHorizon(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);
int IOWriteStatHorizon_dummy(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary);
int IOUpdateProbes(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int IOUpdateProbes_dummy(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int IOWriteProbes(struct DNA_RunOptions *RunOptions);
int IOWriteProbes_dummy(struct DNA_RunOptions *RunOptions);
int IOIdentifySampleIDs(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields);
int IOProgressInitial();
int IOProgressUpdate(int *prog, DNA_FLOAT elapsedtime, DNA_FLOAT totaltime);
int IOProgressFinal();

#endif
