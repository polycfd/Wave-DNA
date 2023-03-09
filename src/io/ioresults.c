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

  RunOptions->WaveResults = fopen("results/fields.dat", "a");
  fprintf(RunOptions->WaveResults, "===============================\n");
  fprintf(RunOptions->WaveResults, "time   ");
  fprintf(RunOptions->WaveResults, "%.6e", RunOptions->t);
  fprintf(RunOptions->WaveResults, "\n");
  fprintf(RunOptions->WaveResults, "timeStep   ");
  fprintf(RunOptions->WaveResults, "%d", RunOptions->NumericsFD.dtNumber);
  fprintf(RunOptions->WaveResults, "\n");
  fprintf(RunOptions->WaveResults, "R   ");
  fprintf(RunOptions->WaveResults, "%.6e", MovingBoundary->R);
  fprintf(RunOptions->WaveResults, "\n");
  fprintf(RunOptions->WaveResults, "dotR   ");
  fprintf(RunOptions->WaveResults, "%.6e", MovingBoundary->U);
  fprintf(RunOptions->WaveResults, "\n");
  for (iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints); iPoint++)
  {
    fprintf(RunOptions->WaveResults, "%d", iPoint);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "%.6e", Fields->Grid.x[iPoint]);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "%.6e", Fields->PressureField.val[iPoint]);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "%.6e", Fields->phi.val[iPoint]);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "%.6e", Fields->BackgroundVelocity.val[iPoint]);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "\n");
  }
  fprintf(RunOptions->WaveResults, "EOS\n");
  fprintf(RunOptions->WaveResults, "===============================\n");
  fclose(RunOptions->WaveResults);

  return 0;
}

int IOWriteProbes(struct DNA_RunOptions *RunOptions)
{
  RunOptions->WaveResults = fopen("results/probes.dat", "a");
  for (int iTime = 0; iTime < (RunOptions->writeFrequency); iTime++)
  {
    fprintf(RunOptions->WaveResults, "%d", RunOptions->Probes.TimeIDs[iTime]);
    fprintf(RunOptions->WaveResults, "   ");
    fprintf(RunOptions->WaveResults, "%.6e", RunOptions->Probes.SampleTime[iTime]);
    fprintf(RunOptions->WaveResults, "   ");
    for (int i = 0; i < (RunOptions->Probes.nSamplePoints); i++)
    {
      fprintf(RunOptions->WaveResults, "%.6e", RunOptions->Probes.SamplePressure[iTime][i]);
      fprintf(RunOptions->WaveResults, "   ");
    }
    fprintf(RunOptions->WaveResults, "\n");
  }
  fclose(RunOptions->WaveResults);

  return 0;
}

int IOWriteProbes_dummy(struct DNA_RunOptions *RunOptions) { return 0; }

int IOWriteHeader(struct DNA_RunOptions *RunOptions)
{
  RunOptions->WaveResults = fopen("results/fields.dat", "a");
  fprintf(RunOptions->WaveResults, "iPoint; x[m]; p1[Pa]; phi1[m^2/s]; u0[m/s]\n\n");
  fclose(RunOptions->WaveResults);

  if (RunOptions->Probes.nSamplePoints > 0)
  {
    RunOptions->WaveResults = fopen("results/probes.dat", "a");
    fprintf(RunOptions->WaveResults, "timeStep; time[s]; ");
    for (int i = 0; i < (RunOptions->Probes.nSamplePoints); i++)
    {
      fprintf(RunOptions->WaveResults, "p1(x=");
      fprintf(RunOptions->WaveResults, "%.6e", RunOptions->Probes.SamplePoints[i]);
      fprintf(RunOptions->WaveResults, "m)[Pa]; ");
    }
    fprintf(RunOptions->WaveResults, "\n");
    fclose(RunOptions->WaveResults);
  }

  if (RunOptions->BoundaryMotionType == 2 || RunOptions->BoundaryMotionType == 3)
  {
    RunOptions->WaveResults = fopen("results/horizon.dat", "a");
    fprintf(RunOptions->WaveResults, "# timeStep; time[s]; R[m]; dotR[m/s]; rh[m]; p1[Pa]; phi1[m^2/s]");
    fprintf(RunOptions->WaveResults, "\n");
    fclose(RunOptions->WaveResults);
  }

  return 0;
}

int IOWriteStatHorizon(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary)
{
  if (RunOptions->statHorizon == 1)
  {
    DNA_FLOAT rl;
    DNA_FLOAT rr;
    DNA_FLOAT rh;
    DNA_FLOAT pl;
    DNA_FLOAT pr;
    DNA_FLOAT prh;
    DNA_FLOAT phil;
    DNA_FLOAT phir;
    DNA_FLOAT phirh;
    DNA_FLOAT horizonIndicator;

    for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints - 1); iPoint++)
    {
      horizonIndicator = (RunOptions->radius_horizon - Fields->Grid.x[iPoint]) * (RunOptions->radius_horizon - Fields->Grid.x[iPoint + 1]);

      if (horizonIndicator < 0.0)
      {
        rh = RunOptions->radius_horizon;
        rl = Fields->Grid.x[iPoint];
        rr = Fields->Grid.x[iPoint + 1];

        pl = Fields->PressureField.val[iPoint];
        pr = Fields->PressureField.val[iPoint + 1];
        prh = (pr - pl) / (rr - rl) * (rh - rl) + pl;

        phil = Fields->phi.val[iPoint];
        phir = Fields->phi.val[iPoint + 1];
        phirh = (phir - phil) / (rr - rl) * (rh - rl) + phil;

        RunOptions->WaveResults = fopen("results/horizon.dat", "a");
        fprintf(RunOptions->WaveResults, "%d", RunOptions->NumericsFD.dtNumber);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", RunOptions->t);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", MovingBoundary->R);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", MovingBoundary->U);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", rh);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", prh);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "%.6e", phirh);
        fprintf(RunOptions->WaveResults, "   ");
        fprintf(RunOptions->WaveResults, "\n");
        fclose(RunOptions->WaveResults);
      }
    }
  }

  return 0;
}

int IOWriteStatHorizon_dummy(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, struct DNA_MovingBoundary *MovingBoundary) { return 0; }

// Thi function identiies the sample IDs based on the user-specific sample locations
int IOIdentifySampleIDs(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields)
{
  int success;
  DNA_FLOAT xl;
  DNA_FLOAT xr;
  DNA_FLOAT dxl;
  DNA_FLOAT dxr;
  DNA_FLOAT xtarg;

  for (int iSample = 0; iSample < (RunOptions->Probes.nSamplePoints); iSample++)
  {
    xtarg = RunOptions->Probes.SamplePoints[iSample];
    success = 0;

    for (int iPoint = 0; iPoint < (RunOptions->NumericsFD.NPoints - 1); iPoint++)
    {
      xl = Fields->Grid.x[iPoint];
      xr = Fields->Grid.x[iPoint + 1];

      if ((xl <= xtarg) & (xr >= xtarg))
      {
        dxl = xtarg - xl;
        dxr = xr - xtarg;

        if (dxl < dxr)
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
      printf("+ WARNING: sample point no %d has no proper nearest neighbour!\n", iSample + 1);

      if (xtarg < Fields->Grid.xmov)
      {
        RunOptions->Probes.SampleIDs[iSample] = 0;
      }
      else
      {
        RunOptions->Probes.SampleIDs[iSample] = RunOptions->NumericsFD.NPoints - 1;
      }
    }
  }

  printf("+ Pressure probes are taken at:\n");
  printf("+ IDs:   ");
  for (int iSample = 0; iSample < (RunOptions->Probes.nSamplePoints); iSample++)
  {
    printf("%d", RunOptions->Probes.SampleIDs[iSample]);
    printf("   ");
  }
  printf("\n");
  printf("+ Target locations:   ");
  for (int iSample = 0; iSample < (RunOptions->Probes.nSamplePoints); iSample++)
  {
    printf("%.6e", RunOptions->Probes.SamplePoints[iSample]);
    printf("   ");
  }
  printf("\n");
  printf("+ Actual locations:   ");
  for (int iSample = 0; iSample < (RunOptions->Probes.nSamplePoints); iSample++)
  {
    printf("%.6e", Fields->Grid.x[RunOptions->Probes.SampleIDs[iSample]]);
    printf("   ");
  }
  printf("\n");

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
    RunOptions->Probes.SamplePressure[id][i] = Fields->PressureField.val[RunOptions->Probes.SampleIDs[i]];
  }

  return 0;
}

int IOUpdateProbes_dummy(int id, struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields) { return 0; }
