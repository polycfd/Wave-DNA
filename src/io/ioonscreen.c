#include "DNA.h"
#include "DNA-functions.h"

int IOInfoOnScreen(int Priority, int PriorityOption, char* str)
{
  if (Priority & PriorityOption) printf("+ %s\n", str);
  return 0;
}

int IOErrorOnScreen(int num, char* message)
{
  char str[DNA_STRINGLENGTH_SPRINTF];

  if (num != 0) /* Terminal error. Stop program and exit */
  {
    printf("+\n");
    printf("+ --------------------------------------------------------------------------- \n");
    sprintf(str, "ERROR: %s", message);
    IOInfoOnScreen(DNA_HIGHPRIO, DNA_HIGHPRIO, str);
    printf("+ --------------------------------------------------------------------------- \n");
    exit(1);
  }
  else
  {
    printf("+ --------------------------------------------------------------------------- \n");
    sprintf(str, "WARNING: %s", message);
    IOInfoOnScreen(DNA_HIGHPRIO, DNA_HIGHPRIO, str);
    printf("+ --------------------------------------------------------------------------- \n+\n");
  }

  return 0;
}

int IOWelcomeScreen()
{
  printf("+ --------------------------------------------------------------------------- \n");
  printf("+ \n");
  printf("+ Name... \n");
  printf("+ \n");
  printf("+ Version: %.1f \n", DNA_VERSION_NUM);
  printf("+ Release date: %s  \n", DNA_RELEASE_DATE);
#ifdef __OPTIMIZE__
  printf("+ Compiled in release mode \n");
#else
  printf("+ Compiled in debug mode \n");
#endif
#if defined(DNA_PRECISION_QUAD)
  printf("+ Compiled with quad precision \n");
#elif defined(DNA_PRECISION_LONGDOUBLE)
  printf("+ Compiled with long-double precision \n");
#else
  printf("+ Compiled with double precision \n");
#endif
  printf("+\n");
  printf("+ --------------------------------------------------------------------------- \n+\n");

  return 0;
}

int IOExitScreen(double totaltime)
{
  printf("+ --------------------------------------------------------------------------- \n");
  printf("+ \n");
  printf("+ Total execution time: %.3f s.\n", totaltime);
  printf("+                                                                    The END! \n");
  printf("+ --------------------------------------------------------------------------- \n");

  return 0;
}
