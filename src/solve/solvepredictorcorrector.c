#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions are related to the predictor-corrector method. In order to
compute the predicted and corrected solutions for the acoustic perturbation potential
phi, the functions make use of the fields AA, BB, and RHS, which are computed in
<transformeqn.c>.
---------------------------------------------------------**/

/** Re-evaluate the Laplacian term needed for the corrected solution. The gradient
is required as the Laplacian for a wave with variable surface area of the (curved)
wavefront involves a gradient contribution as well. **/
int SolveCorrectLaplacian(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  FDGradient(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI1_phi);
  FDLaplacian(&RunOptions->NumericsFD, &Fields->phi, &Fields->dXI2_phi);

  return 0;
}

/** Predicted solution == final solution if the number of corrector steps is zero. **/
int SolveExplicit_Predictor(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;

  for (iPoint = 1; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->phi.val[iPoint] = (Fields->RHS.val[iPoint] + Fields->BB.val[iPoint]) / Fields->AA.val[iPoint];
  }

  return 0;
}

/** Corrected solution based on the re-evaluated Laplacian. The corrected solution
consists of the initial guess and the solution obtained for the field based on the
re-evaluated Laplacian, where the two contributions are weighted by (1-gamma) and
gamma, respectively, with gamma being the corrector weight as specified in the 
options file. Typically 0 <= gamma <= 1, but gamma > 1 is admissible as well. **/
int SolveExplicit_Corrector(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;
  DNA_FLOAT gamma = RunOptions->NumericsFD.correctorWeight;

  for (iPoint = 1; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->phi.val[iPoint] =
        (1.0 - gamma) * Fields->phi1_initGuess.val[iPoint] + gamma * Fields->BBbyAA.val[iPoint] + gamma * Fields->RHS.val[iPoint] / Fields->AA.val[iPoint];
  }

  return 0;
}

/** The initial guess, more specifically the term BBbyAA (BB/AA) needs to be stored
in order to compute the corrected solution based on the weighted contributions
constituted by the initial guess and the re-evaluated field based on the re-evaluated
Laplacian. **/
int SolveStoreInitialGuess(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int iPoint;

  for (iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    Fields->BBbyAA.val[iPoint] = Fields->BB.val[iPoint] / (Fields->AA.val[iPoint]);
    Fields->phi1_initGuess.val[iPoint] = Fields->phi.val[iPoint];
  }

  return 0;
}
