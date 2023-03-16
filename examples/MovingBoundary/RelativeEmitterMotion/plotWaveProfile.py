import numpy as np
import matplotlib.pyplot as plt
import sys

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

optionsfile = 'run.DNA' # default if no command line option (plot *.py optionsfile) is given
if len(sys.argv) > 1:
    optionsfile = sys.argv[1]
fo = open(optionsfile, "r")
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
fo.close()

def _read(filename, timestep):

    fo = open(filename, 'r')
    x = []
    p = []

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
                    line = next(fo)
    fo.close()
    return x,p

filename = "results/fields.dat"

x,p = _read(filename, 48000)
plt.plot(x,p, linewidth=4.0, color='red', label='$t=3/f_{\mathrm{a}}$')
x,p = _read(filename, 144000)
plt.plot(x,p, linewidth=2.0, color='black', label='$t=9f_{\mathrm{a}}$')

plt.xlabel('$x/\lambda_0$', fontsize=fontsize)
plt.ylabel('$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)
plt.legend(loc='upper center', frameon=False, ncol=2, fontsize=fontsize, bbox_to_anchor=[0.5,1.2])

#plt.xlim(0,12)
#plt.ylim(-1.1,1.1)

plt.grid()
plt.show()
#plt.savefig("../../figures/RelativeEmitterMotion.pdf", bbox_inches = "tight")
#plt.close()

