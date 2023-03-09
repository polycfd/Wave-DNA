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

  printf("+ Reading options file %s\n+\n", RunOptions.OptionsFile);
  IOReadOptionsFile(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  /** Initialize the wave  **/
  InitializeSimulation(&RunOptions, &Fields, &MovingBoundary);
  InitializeProcessOptions(&RunOptions, &Fields, &MovingBoundary);

  /** Time loop of the numerical simulation **/
  printf("+ Executing solver.\n+\n");
  SolveTimeLoop(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);
  printf("+ Solver has concluded in %.3f s.\n+\n", (double) (clock() - starttime) / CLOCKS_PER_SEC);

  /** Free all allocated memory **/
  MemoryFreeFields(&RunOptions, &Fields);

  IOExitScreen((double) (clock() - starttime) / CLOCKS_PER_SEC);

  return 0;
}
