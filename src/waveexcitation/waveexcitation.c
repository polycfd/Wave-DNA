#include "DNA.h"
#include "DNA-functions.h"

/**--------------------------------------------------------- 
The following functions are meant to impose a wave exciation signal for the
acoustic pressure. The return value (acoustic pressure) is
applied to the excitation node in in the time loop (function <SolveTimeLoop> 
in <solveloop.c>).

The excitation function is either sinusoidal (<WaveExcitationSineModulated>)
or constant (<WaveExcitationConst>). In addition, it can be convoluted with
a Gaussian envelope (<WaveExcitationGaussEnvelope> or
<WaveExcitationGaussEnvelopeNone> if no envelope is applied). The corresponding
excitation and Gauss envelope functions are called by reference based on the 
user settings in the options file. The respective function pointers are called
<PeriodicPressureExcitation> and <GaussConvolutionFactor>.

The function pointers are called in the time loop function <SolveTimeLoop>
(see <solveloop.c>). The sequence in <SolveTimeLoop> is as follows:

1.) Determine the Gauss envelope (= 1.0 = const if no envelope is applied)
2.) Determine the excitation pressure <pExcitation>.

The excitation pressure is then converted into the excitation potential and
assigned to the excitation node in the <Solve> function called inside
<SolveTimeLoop>. The conversion and assignment is done by the function
<TransformExcitationPressureToExcitationPotential>. 
---------------------------------------------------------**/

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
