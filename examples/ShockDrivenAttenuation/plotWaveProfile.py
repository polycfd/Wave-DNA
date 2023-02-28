import numpy as np
import matplotlib.pyplot as plt
import math


#===================================================================
preamble='\\usepackage{amsmath}\n\\usepackage{amsfonts}'

#plt.rc('text',usetex=True)
#plt.rc('text.latex',preamble=preamble)
plt.rcParams['mathtext.fontset']='cm' 

fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6,3))
useMathText=True
fontsize = 14
font = {'family' : 'serif', 'size' : fontsize}
ax.xaxis.label.set_size(fontsize)
ax.yaxis.label.set_size(fontsize)
plt.setp(ax.get_xticklabels(), fontsize=fontsize)
plt.setp(ax.get_yticklabels(), fontsize=fontsize)
#===================================================================

fo = open("run.DNA", "r")
for line in fo:
  if "ExcitationFrequency" in line:
    splitline = line.split()
    fa = float(splitline[1])
  if "PressureAmplitude" in line:
    splitline = line.split()
    pa = float(splitline[1])
  if "SoundSpeed" in line:
    splitline = line.split()
    c0 = float(splitline[1])
  if "Density" in line:
    splitline = line.split()
    rho0 = float(splitline[1])
  if "FluidNonLinearity" in line:
    splitline = line.split()
    beta = float(splitline[1])
fo.close()

lambdaa = c0/fa
xshbylambda = rho0*c0**3/(2.0*math.pi*beta*fa*pa)/lambdaa

def _read(filename, timestep):

    fo = open(filename, 'r')
    x = []
    p = []
    fay = []

    for line in fo:

        if 'timeStep' in line:
            splitline = line.split()
            timeStep = int(splitline[1])

            if timeStep == timestep:
                line = next(fo)
                line = next(fo)
                line = next(fo)

                while "EOS" not in line:
                    splitline = line.split()
                    x.append(float(splitline[1])*fa/c0)
                    p.append(float(splitline[2])/pa)
                    fay.append(math.pi/(1.0 + x[-1]/xshbylambda))
                    line = next(fo)
    fo.close()
    return x,p,fay

filename = "results/fields.dat"

x,p,fay = _read(filename, 384000)
plt.plot(x,p, linewidth=2.0, color='black', linestyle='solid', label='$\mathrm{wave\:profile}$')
plt.plot(x,fay, linewidth=2.0, color='red', linestyle='dashed', label='$\mathrm{Fay\:envelope}$')

plt.xlabel('$x/\lambda_0$', fontsize=fontsize)
plt.ylabel('$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)
plt.legend(loc='lower right', frameon=True, ncol=1, fontsize=fontsize)

plt.xlim(-5,55)
#plt.xlim(30,40)
plt.ylim(-1.5,1.5)

plt.grid()
plt.show()
#plt.savefig("../figures/shock_attenuation.pdf", bbox_inches = "tight")
#plt.savefig("../figures/shock_attenuation_detail.pdf", bbox_inches = "tight")
#plt.close()

