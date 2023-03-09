#include "DNA.h"
#include "DNA-functions.h"

int IOReadCommandLineOptions(int argc, char **args, struct DNA_RunOptions *RunOptions)
{
  char str[DNA_STRINGLENGTH_SPRINTF];

  sprintf(RunOptions->OptionsFile, "./run.DNA");

  int i = 1;  // First argument is the DNA call
  while (i < argc)
  {
    if (strcmp("-options", args[i]) == 0)
    {
      sprintf(RunOptions->OptionsFile, "%s", args[i + 1]);
      i += 2;
    }
    else
    {
      sprintf(str, "Unknown command line options: %s", args[i]);
      IOErrorOnScreen(1, str);
      ++i;
    }
  }
  return 0;
}
