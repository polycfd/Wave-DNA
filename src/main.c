#include "DNA.h"
#include "DNA-functions.h"

int main(int argc, char **args)
{
  /** See <DNA.h> for the declaration and content of the following structure.
  They may contain nested structures. **/
  struct DNA_RunOptions RunOptions;
  struct DNA_Fields Fields;
  struct DNA_MovingBoundary MovingBoundary;
  struct DNA_FluidProperties FluidProperties;

  IOWelcomeScreen();

  clock_t starttime = clock();

  /** Read command line options **/
  IOReadCommandLineOptions(argc, args, &RunOptions);

  /** Run and process options **/
  IODefaultOptions(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  printf("+ Reading options file %s\n+\n", RunOptions.OptionsDir);

  IOReadOptionsFile(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  clock_t startTime = clock();

  /** Initialize the wave  **/
  InitializeSimulation(&RunOptions, &Fields, &MovingBoundary);
  InitializeProcessOptions(&RunOptions, &Fields, &MovingBoundary);

  printf("+ Executing solver.\n+\n");

  /** Time loop of the numerical simulation **/
  SolveTimeLoop(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  /** Free all allocated memory **/
  MemoryFreeFields(&RunOptions, &Fields);

  printf("+ Solver has concluded in %.3f s.\n+\n", (double) (clock() - startTime) / CLOCKS_PER_SEC);

  IOExitScreen((double) (clock() - starttime) / CLOCKS_PER_SEC);

  return 0;
}
