#include "strings.h"
#include <sys/stat.h>
#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
Reading the options file. In case of doubt, this file can always be consulted
to check for valid options.
---------------------------------------------------------**/

int IOReadOptionsFile(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary,
                      struct DNA_FluidProperties *FluidProperties)
{
  int l = 0;
  int line = 0;
  FILE *OptionsFile;
  char str[DNA_STRINGLENGTH_SPRINTF];
  char option[DNA_STRINGLENGTH], option2[DNA_STRINGLENGTH], option3[DNA_STRINGLENGTH];
  int ReadingStatus = 1;
  int KeywordStatus = 1;
  char *dummy;  // Dummy required for string-to-float macro

  if ((OptionsFile = fopen(RunOptions->OptionsDir, "r")) == (FILE *) NULL)
  {
    sprintf(str, "File %s cannot be opened for reading.\n", RunOptions->OptionsDir);
    IOErrorOnScreen(1, str);
  }

  while ((l = IOReadOneOption(OptionsFile, option)) != EOF && ReadingStatus == 1)
  {
    line += l;
    /*****************************************************************************/
    /* KEYWORD RUNOPTIONS                                                               */
    /*****************************************************************************/
    if (strncasecmp(option, "RUNOPTIONS", 10) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "WriteFrequency", 14) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->writeFrequency = atoi(option3);
        }
        else if (strncasecmp(option2, "SamplePoints", 12) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->Probes.nSamplePoints = atoi(option3);

          if (RunOptions->Probes.nSamplePoints > 0)
          {
            MemoryAllocSamplePoints(RunOptions->Probes.nSamplePoints, &RunOptions->Probes);

            for (int i = 0; i < RunOptions->Probes.nSamplePoints; i++)
            {
              l = IOReadOneOption(OptionsFile, option3);
              RunOptions->Probes.SamplePoints[i] = DNA_STRINGTOFLOAT(option3);
            }
          }
        }
        else if (strncasecmp(option2, "timestart", 9) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->tStart = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "timeend", 7) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->tEnd = DNA_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of RUNOPTIONS is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /* KEYWORD FLUID                                                         */
    /*****************************************************************************/
    else if (strncasecmp(option, "FLUID", 5) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "SoundSpeed", 10) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          FluidProperties->c0 = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "Density", 7) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          FluidProperties->rho0 = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "FluidNonLinearity", 17) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          FluidProperties->beta = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "LocalNonlinearity", 17) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "True", 4) == 0)
            {
              RunOptions->LagrangianDensityFlag = 1.0;
            }
            else if (strncasecmp(option3, "False", 5) == 0)
            {
              RunOptions->LagrangianDensityFlag = 0.0;
            }
          }
        }
        else
        {
          sprintf(str, "An unknown option of FLUID is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /* KEYWORD EXCITATION                                                         */
    /*****************************************************************************/
    else if (strncasecmp(option, "EXCITATION", 10) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "ExcitationNode", 14) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->iExcitation = atoi(option3);
        }
        else if (strncasecmp(option2, "GaussEnvelope", 13) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "False", 5) == 0)
            {
              RunOptions->WaveExcitation.GaussEnvelopeOption = 0;
            }
            else if (strncasecmp(option3, "True", 4) == 0)
            {
              RunOptions->WaveExcitation.GaussEnvelopeOption = 1;
            }
          }
        }
        else if (strncasecmp(option2, "ExcitationFunctionType", 22) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "sine", 4) == 0)
            {
              RunOptions->WaveExcitation.ExcitationFunction = 0;
            }
            else if (strncasecmp(option3, "GaussPulse", 10) == 0)
            {
              RunOptions->WaveExcitation.ExcitationFunction = 1;
            }
          }
        }
        else if (strncasecmp(option2, "ExcitationStart", 15) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.ExcitationStart = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "ExcitationEnd", 13) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.ExcitationEnd = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "PressureAmplitude", 17) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.PressureAmplitude = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "ExcitationFrequency", 19) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.ExcitationFrequency = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "PhaseShiftAngle", 15) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.PhaseShiftAngle = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "PowerCoeff", 10) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.PowerCoeff = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "EnvelopedPeriod", 15) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.EnvelopedPeriod = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "GaussShapeCoeff", 15) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.GaussShapeCoeff = DNA_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of EXCITATION is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /* KEYWORD BOUNDARYMOTION                                                         */
    /*****************************************************************************/
    else if (strncasecmp(option, "BOUNDARYMOTION", 14) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "BoundaryMotionType", 18) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "linear", 6) == 0)
            {
              RunOptions->BoundaryMotionType = 0;
            }
            else if (strncasecmp(option3, "oscillating", 11) == 0)
            {
              RunOptions->BoundaryMotionType = 1;
            }
            else if (strncasecmp(option3, "stationaryBlackHole", 19) == 0)
            {
              RunOptions->BoundaryMotionType = 2;
            }
            else if (strncasecmp(option3, "stationaryWhiteHole", 19) == 0)
            {
              RunOptions->BoundaryMotionType = 3;
            }
          }
        }
        else if (strncasecmp(option2, "BoundaryMotionStartTime", 23) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.BoundaryMotionStartTime = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "BoundaryMotionEndTime", 21) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.BoundaryMotionEndTime = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "InitialMovingBoundaryPosition", 27) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          MovingBoundary->R0 = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "FixedBoundaryPosition", 21) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          Fields->Grid.xfix = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "MovingBoundaryVelocityAmplitude", 29) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.MovingBoundaryVelocityAmplitude = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "MovingBoundaryFrequency", 21) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->WaveExcitation.MovingBoundaryFrequency = DNA_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of BOUNDARYMOTION is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /* KEYWORD BACKGROUNDFLOW                                                         */
    /*****************************************************************************/
    else if (strncasecmp(option, "BACKGROUNDFLOW", 14) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "BackgroundMotionMode", 20) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "quiescent", 9) == 0)
            {
              RunOptions->backgroundMotionMode = 0;
            }
            else if (strncasecmp(option3, "coupledToMovingBoundary", 23) == 0)
            {
              RunOptions->backgroundMotionMode = 1;
            }
          }
        }
        else if (strncasecmp(option2, "Geometry", 8) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "1dCartesian", 11) == 0)
            {
              RunOptions->dimension = 1;
            }
            else if (strncasecmp(option3, "3dSphericallySymmetric", 22) == 0)
            {
              RunOptions->dimension = 3;
            }
          }
        }
        else if (strncasecmp(option2, "HorizonRadius", 13) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->radius_horizon = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "BackgroundVelocityScalingFactor", 31) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->U_backgroundScaling = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "GravitationalPotential", 22) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "False", 5) == 0)
            {
              RunOptions->GravitationalPotentialOption = 0;
            }
            else if (strncasecmp(option3, "True", 4) == 0)
            {
              RunOptions->GravitationalPotentialOption = 1;
            }
          }
        }
        else
        {
          sprintf(str, "An unknown option of BACKGROUNDFLOW is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /* KEYWORD FINITEDIFFERENCE                                                         */
    /*****************************************************************************/
    else if (strncasecmp(option, "FINITEDIFFERENCE", 16) == 0)
    {
      KeywordStatus = 1;
      while (KeywordStatus == 1 && (l = IOReadOneOption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          KeywordStatus = 0;
        }
        else if (strncasecmp(option2, "NPoints", 7) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->NumericsFD.NPoints = atoi(option3);
        }
        else if (strncasecmp(option2, "TimeStepSize", 12) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->NumericsFD.dt = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "BoundaryConditionEast", 21) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "scattering", 10) == 0)
            {
              RunOptions->NumericsFD.BCEast = 0;
            }
            else if (strncasecmp(option3, "absorbing", 9) == 0)
            {
              RunOptions->NumericsFD.BCEast = 1;
            }
          }
        }
        else if (strncasecmp(option2, "BoundaryConditionWest", 21) == 0)
        {
          if ((l = IOReadOneOption(OptionsFile, option3)) == EOF)
          {
            KeywordStatus = 0;
            ReadingStatus = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "scattering", 10) == 0)
            {
              RunOptions->NumericsFD.BCWest = 0;
            }
            else if (strncasecmp(option3, "absorbing", 9) == 0)
            {
              RunOptions->NumericsFD.BCWest = 1;
            }
          }
        }
        else if (strncasecmp(option2, "correctorWeight", 15) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->NumericsFD.correctorWeight = DNA_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "nCorrectors", 11) == 0)
        {
          l = IOReadOneOption(OptionsFile, option3);
          RunOptions->NumericsFD.nCorrectors = atoi(option3);
        }
        else
        {
          sprintf(str, "An unknown option of FINITEDIFFERENCE is given: %s, line %i", option2, line);
          IOErrorOnScreen(1, str);
        }
      }
    }
    /*****************************************************************************/
    /*                                                                           */
    /*****************************************************************************/
    else
    {
      sprintf(str, "An unknown Keyword given: %s, line %i", option, line);
      IOErrorOnScreen(1, str);
    }
  }
  /* Continue Reading File */
  /* End of reading options file */
  fclose(OptionsFile);

  return 0;
}
