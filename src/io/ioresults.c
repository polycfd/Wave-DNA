#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions are called to write the output files

fields.dat (always created)
probes.dat (optional)
horizon.dat (optional)

fields.dat: Containts the field data such as acoustic pressure over the spatial
            coordinate. The output is written at the write frequency as specified
            in the options file.

probes.dat: Acoustic pressure over time recorded at the grid IDs that correspond
            to the locations specified by the users (nearest neighbors (NB) are
            identified). Note that the corresponding grid points may be moving.
            In the current implementation, no interpolation is done. The NB Ids
            associated with a certain fixed location would likely change in the
            case of considerable grid deformation, which would requires a 
            computationally efficient function to avoid the call of a costly 
            search algorithm at each write instance, in particular if multiple
            sample points are specified.

horizon.dat:  Time signals taken at the sonic horizon of a steady-state acoustic
              black or white hole. A linear interpolation between the neighbor
              grid points of the horizon is carried out. Here, we have to apply 
              a reltively costly search algorithm to identify the neighbors, which
              typically change over time due to the dilation of the physical
              domain. The function could be improved in terms of speed.
---------------------------------------------------------**/

int IOWriteToDisc(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{
  int iPoint;

  RunOptions->Results = fopen("results/fields.dat", "a");
  fprintf(RunOptions->Results, "===============================\n");
  fprintf(RunOptions->Results, "time   ");
  fprintf(RunOptions->Results, "%.6e", RunOptions->t);
  fprintf(RunOptions->Results, "\n");
  fprintf(RunOptions->Results, "timeStep   ");
  fprintf(RunOptions->Results, "%d", RunOptions->NumericsFD.dtNumber);
  fprintf(RunOptions->Results, "\n");
  fprintf(RunOptions->Results, "R   ");
  fprintf(RunOptions->Results, "%.6e", MovingBoundary->R);
  fprintf(RunOptions->Results, "\n");
  fprintf(RunOptions->Results, "dotR   ");
  fprintf(RunOptions->Results, "%.6e", MovingBoundary->U);
  fprintf(RunOptions->Results, "\n");
  for (iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    fprintf(RunOptions->Results, "%d", iPoint);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "%.6e", Fields->Grid.x[iPoint]);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "%.6e", Fields->PhiField.PressureField[iPoint]);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "%.6e", Fields->PhiField.phi[iPoint]);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "%.6e", Fields->BackgroundFlowField.BackgroundVelocity[iPoint]);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "\n");
  }
  fprintf(RunOptions->Results, "EOS\n");
  fprintf(RunOptions->Results, "===============================\n");
  fclose(RunOptions->Results);

  return 0;
}

int IOWriteProbes(struct DNA_RunOptions *RunOptions)
{
  RunOptions->Results = fopen("results/probes.dat", "a");
  for (int iTime = 0; iTime < (RunOptions->writeFrequency); iTime++)
  {
    fprintf(RunOptions->Results, "%d", RunOptions->Probes.TimeIDs[iTime]);
    fprintf(RunOptions->Results, "   ");
    fprintf(RunOptions->Results, "%.6e", RunOptions->Probes.SampleTime[iTime]);
    fprintf(RunOptions->Results, "   ");
    for (int i = 0; i < (RunOptions->Probes.nSamplePoints); i++)
    {
      fprintf(RunOptions->Results, "%.6e", RunOptions->Probes.SamplePressure[iTime][i]);
      fprintf(RunOptions->Results, "   ");
    }
    fprintf(RunOptions->Results, "\n");
  }
  fclose(RunOptions->Results);

  return 0;
}

int IOWriteProbes_dummy(struct DNA_RunOptions *RunOptions) { return 0; }

int IOWriteHeader(struct DNA_RunOptions *RunOptions)
{
  RunOptions->Results = fopen("results/fields.dat", "a");
  fprintf(RunOptions->Results, "iPoint; x[m]; p1[Pa]; phi1[m^2/s]; u0[m/s]\n\n");
  fclose(RunOptions->Results);

  if (RunOptions->Probes.nSamplePoints > 0)
  {
    RunOptions->Results = fopen("results/probes.dat", "a");
    fprintf(RunOptions->Results, "timeStep; time[s]; ");
    for (int i = 0; i < (RunOptions->Probes.nSamplePoints); i++)
    {
      fprintf(RunOptions->Results, "p1(x=");
      fprintf(RunOptions->Results, "%.6e", RunOptions->Probes.SamplePoints[i]);
      fprintf(RunOptions->Results, "m)[Pa]; ");
    }
    fprintf(RunOptions->Results, "\n");
    fclose(RunOptions->Results);
  }

  if (RunOptions->BoundaryMotionType == 2 || RunOptions->BoundaryMotionType == 3)
  {
    RunOptions->Results = fopen("results/horizon.dat", "a");
    fprintf(RunOptions->Results, "# timeStep; time[s]; R[m]; dotR[m/s]; rh[m]; p1[Pa]; phi1[m^2/s]");
    fprintf(RunOptions->Results, "\n");
    fclose(RunOptions->Results);
  }

  return 0;
}

int IOWriteStatHorizon(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{
  if (RunOptions->statHorizon == 1)
  {
    for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints - 1); iPoint++)
    {
      DNA_FLOAT horizonIndicator = (RunOptions->radius_horizon - Fields->Grid.x[iPoint]) * (RunOptions->radius_horizon - Fields->Grid.x[iPoint + 1]);

      if (horizonIndicator < 0.0)
      {
        DNA_FLOAT rh = RunOptions->radius_horizon;
        DNA_FLOAT rl = Fields->Grid.x[iPoint];
        DNA_FLOAT rr = Fields->Grid.x[iPoint + 1];

        DNA_FLOAT pl = Fields->PhiField.PressureField[iPoint];
        DNA_FLOAT pr = Fields->PhiField.PressureField[iPoint + 1];
        DNA_FLOAT prh = (pr - pl) / (rr - rl) * (rh - rl) + pl;

        DNA_FLOAT phil = Fields->PhiField.phi[iPoint];
        DNA_FLOAT phir = Fields->PhiField.phi[iPoint + 1];
        DNA_FLOAT phirh = (phir - phil) / (rr - rl) * (rh - rl) + phil;

        RunOptions->Results = fopen("results/horizon.dat", "a");
        fprintf(RunOptions->Results, "%d", RunOptions->NumericsFD.dtNumber);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", RunOptions->t);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", MovingBoundary->R);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", MovingBoundary->U);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", rh);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", prh);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "%.6e", phirh);
        fprintf(RunOptions->Results, "   ");
        fprintf(RunOptions->Results, "\n");
        fclose(RunOptions->Results);
      }
    }
  }

  return 0;
}

int IOWriteStatHorizon_dummy(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary) { return 0; }

// Thi function identiies the sample IDs based on the user-specific sample locations
int IOIdentifySampleIDs(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  if (RunOptions->Probes.nSamplePoints)
  {
    for (int iSample = 0; iSample < RunOptions->Probes.nSamplePoints; iSample++)
    {
      DNA_FLOAT xtarg = RunOptions->Probes.SamplePoints[iSample];
      int success = 0;

      for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints - 1); iPoint++)
      {
        DNA_FLOAT x1 = Fields->Grid.x[iPoint];
        DNA_FLOAT x2 = Fields->Grid.x[iPoint + 1];

        DNA_FLOAT dx1 = x1 - xtarg;
        DNA_FLOAT dx2 = x2 - xtarg;

        // If xtarg is betwen x1 and x2, then dx1*dx2 < 0.
        // By checking for dx1*dx2 < DNA_SMALL, boundary points are identified as well (xtarg = x1 or xtarg = x2).
        if (dx1 * dx2 < DNA_SMALL)
        {
          if (DNA_ABS(dx1) < DNA_ABS(dx2))
          {
            RunOptions->Probes.SampleIDs[iSample] = iPoint;
          }
          else
          {
            RunOptions->Probes.SampleIDs[iSample] = iPoint + 1;
          }
          success = 1;
        }
      }

      if (!success)
      {
        printf("+ INFO: sample point no %d is a domain boundary\n", iSample + 1);

        DNA_FLOAT xmin = APECSS_MIN(Fields->Grid.xmov, Fields->Grid.xfix);

        int xminID;
        int xmaxID;

        if (Fields->Grid.x[RunOptions->NumericsFD.NPoints - 1] > Fields->Grid.x[0])
        {
          xminID = 0;
          xmaxID = RunOptions->NumericsFD.NPoints - 1;
        }
        else
        {
          xminID = RunOptions->NumericsFD.NPoints - 1;
          xmaxID = 0;
        }

        if (xtarg < xmin)
        {
          RunOptions->Probes.SampleIDs[iSample] = xminID;
        }
        else
        {
          RunOptions->Probes.SampleIDs[iSample] = xmaxID;
        }
      }
    }

    printf("+ Pressure probes are taken at:\n");
    printf("+ IDs:   ");
    for (int iSample = 0; iSample < RunOptions->Probes.nSamplePoints; iSample++)
    {
      printf("%d", RunOptions->Probes.SampleIDs[iSample]);
      printf("   ");
    }
    printf("\n");

    printf("+ Target locations:   ");
    for (int iSample = 0; iSample < RunOptions->Probes.nSamplePoints; iSample++)
    {
      printf("%.6e", RunOptions->Probes.SamplePoints[iSample]);
      printf("   ");
    }
    printf("\n");
    printf("+ Actual locations:   ");
    for (int iSample = 0; iSample < RunOptions->Probes.nSamplePoints; iSample++)
    {
      printf("%.6e", Fields->Grid.x[RunOptions->Probes.SampleIDs[iSample]]);
      printf("   ");
    }
    printf("\n");
  }

  return 0;
}

// This function updates the probe vectors to be written
// to disk at the next write instance
int IOUpdateProbes(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  RunOptions->Probes.SampleTime[id] = RunOptions->t;
  RunOptions->Probes.TimeIDs[id] = RunOptions->NumericsFD.dtNumber;

  for (int i = 0; i < RunOptions->Probes.nSamplePoints; i++)
  {
    RunOptions->Probes.SamplePressure[id][i] = Fields->PhiField.PressureField[RunOptions->Probes.SampleIDs[i]];
  }

  return 0;
}

int IOUpdateProbes_dummy(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields) { return 0; }
