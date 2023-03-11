#include "DNA.h"
#include "DNA-functions.h"

int IOErrorOnScreen(int num, char* message)
{
  if (num != 0) /* Terminal error. Stop program and exit */
  {
    printf("+\n");
    printf("+ --------------------------------------------------------------------------- \n");
    printf("+ ERROR: %s", message);
    printf("+ --------------------------------------------------------------------------- \n");
    exit(1);
  }
  else
  {
    printf("+ --------------------------------------------------------------------------- \n");
    printf("+ WARNING: %s", message);
    printf("+ --------------------------------------------------------------------------- \n+\n");
  }

  return 0;
}

int IOWelcomeScreen()
{
  printf("+ --------------------------------------------------------------------------- \n");
  printf("+ \n");
  printf("+ ██╗    ██╗ █████╗ ██╗   ██╗███████╗    ██████╗ ███╗   ██╗ █████╗ \n");
  printf("+ ██║    ██║██╔══██╗██║   ██║██╔════╝    ██╔══██╗████╗  ██║██╔══██╗ \n");
  printf("+ ██║ █╗ ██║███████║██║   ██║█████╗█████╗██║  ██║██╔██╗ ██║███████║ \n");
  printf("+ ██║███╗██║██╔══██║╚██╗ ██╔╝██╔══╝╚════╝██║  ██║██║╚██╗██║██╔══██║ \n");
  printf("+ ╚███╔███╔╝██║  ██║ ╚████╔╝ ███████╗    ██████╔╝██║ ╚████║██║  ██║ \n");
  printf("+  ╚══╝╚══╝ ╚═╝  ╚═╝  ╚═══╝  ╚══════╝    ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝ \n");
  printf("+ \n");
  printf("+ Version: %.1f \n", DNA_VERSION_NUM);
  printf("+ Release date: %s  \n", DNA_RELEASE_DATE);
#ifdef __OPTIMIZE__
  printf("+ Compiled in release mode \n");
#else
  printf("+ Compiled in debug mode \n");
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

int IOProgressInitial()
{
  fprintf(stdout, "+ Progress %%: ");
  fflush(stdout);
  return (0);
}

int IOProgressUpdate(int* prog, DNA_FLOAT elapsedtime, DNA_FLOAT totaltime)
{
  if (elapsedtime > (DNA_FLOAT) *prog * 0.02 * totaltime)
  {
    if (!(*prog % 5) && *prog != 50)
    {
      fprintf(stdout, "%i", *prog * 2);
      fflush(stdout);
    }
    else
    {
      fprintf(stdout, ".");
      fflush(stdout);
    }

    (*prog)++;
  }

  return (0);
}

int IOProgressFinal()
{
  fprintf(stdout, "100\n+\n");
  return (0);
}