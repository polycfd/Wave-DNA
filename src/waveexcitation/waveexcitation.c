#include "DNA.h"
#include "DNA-functions.h"


int WaveExcitationSetPressure(struct DNA_RunOptions *RunOptions, struct DNA_Fields *Fields, DNA_FLOAT time)
{
  if((time >= RunOptions->WaveExcitation.ExcitationStart) && (time <= RunOptions->WaveExcitation.ExcitationEnd))
  {
    Fields->phi.val[RunOptions->iExcitation] = RunOptions->pExcitation;
  }
  else{
    Fields->phi.val[RunOptions->iExcitation] = 0.0;
  }
  
  return 0;
}

DNA_FLOAT WaveExcitationSineModulated(struct DNA_RunOptions *RunOptions)
{ 
  DNA_FLOAT yamp = RunOptions->WaveExcitation.PressureAmplitude;
  DNA_FLOAT frequ = RunOptions->WaveExcitation.ExcitationFrequency;
  DNA_FLOAT deltaPhi = RunOptions->WaveExcitation.PhaseShiftAngle;
  DNA_FLOAT power = RunOptions->WaveExcitation.PowerCoeff;
  DNA_FLOAT time = RunOptions->t;
  DNA_FLOAT GaussEnvelope = RunOptions->WaveExcitation.GaussConvolutionFactor;
  
  return (yamp*DNA_POW(DNA_SIN((DNA_FLOAT)2.0*DNA_PI*frequ*time - deltaPhi), power))*GaussEnvelope;
}

DNA_FLOAT WaveExcitationConst(struct DNA_RunOptions *RunOptions)
{
  DNA_FLOAT GaussEnvelope = RunOptions->WaveExcitation.GaussConvolutionFactor;
  return RunOptions->WaveExcitation.PressureAmplitude*GaussEnvelope;
}

DNA_FLOAT WaveExcitationGaussEnvelope(DNA_FLOAT periodNo, DNA_FLOAT frequ, DNA_FLOAT time, DNA_FLOAT shapeCoeff)
{
  DNA_FLOAT GaussEnvelope;
  DNA_FLOAT arg = 4.0*DNA_POW2(frequ*time - periodNo + (DNA_FLOAT)0.5);
  GaussEnvelope = DNA_POW(shapeCoeff, arg);
  
  return GaussEnvelope;
}

DNA_FLOAT WaveExcitationGaussEnvelopeNone(DNA_FLOAT periodNo, DNA_FLOAT frequ, DNA_FLOAT time, DNA_FLOAT shapeCoeff)
{ 
  return 1.0;
}
