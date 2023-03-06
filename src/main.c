#include "DNA.h"
#include "DNA-functions.h"

int main(int argc, char **args)
{
  char str[DNA_STRINGLENGTH_SPRINTF];

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

  snprintf(str, sizeof(str), "Reading options file %s\n+", RunOptions.OptionsDir);
  IOInfoOnScreen(DNA_HIGHPRIO, RunOptions.OnScreenIOPriority, str);
  IOReadOptionsFile(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  clock_t startTime = clock();

  IOInfoOnScreen(DNA_HIGHPRIO, RunOptions.OnScreenIOPriority, str);

  /** Initialize the wave  **/
  InitializeSimulation(&RunOptions, &Fields, &MovingBoundary);
  InitializeProcessOptions(&RunOptions, &Fields, &MovingBoundary);

  sprintf(str, "Executing solver.");
  IOInfoOnScreen(DNA_MEDIUMPRIO, RunOptions.OnScreenIOPriority, str);

  /** Time loop of the numerical simulation **/  
  SolveTimeLoop(&RunOptions, &Fields, &MovingBoundary, &FluidProperties);

  /** Free all allocated memory **/
  MemoryFreeFields(&RunOptions, &Fields);

  sprintf(str, "Solver has concluded in %.3f s.\n+", (double) (clock() - startTime) / CLOCKS_PER_SEC);
  IOInfoOnScreen(DNA_HIGHPRIO, RunOptions.OnScreenIOPriority, str);

  IOExitScreen((double) (clock() - starttime) / CLOCKS_PER_SEC);

  return 0;
}
