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

x,p = _read(filename, 200000)
plt.plot(x,p, linewidth=2.0, color='black', label=r'$\mathrm{simulation}$')

decay_pos = []
decay_neg = []
for i in range(0,len(x)):
    decay_pos.append(1.0/x[i])
    decay_neg.append(-1.0/x[i])
plt.plot(x,decay_pos, linewidth=2.0, color='red', linestyle='dashed', label=r'$1/r\mathrm{\:decay}$')
plt.plot(x,decay_neg, linewidth=2.0, color='red', linestyle='dashed')

plt.xlabel(r'$r/\lambda_\mathrm{a}$', fontsize=fontsize)
plt.ylabel(r'$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)

plt.xlim((0,12))
plt.ylim((-1, 1))
plt.xticks([0, 2, 4, 6, 8, 10, 12], ['$0$', '$2$', '$4$', '$6$', '$8$', '$10$', '$12$'], fontsize=fontsize)
plt.yticks([-1, 0, 1], ['$-1$', '$0$', '$1$'], fontsize=fontsize)

plt.legend(fontsize=fontsize, loc='upper right', ncol=1, frameon=True)
plt.grid()

plt.show()
#plt.savefig("wave.pdf", bbox_inches = "tight", dpi=400)
#plt.close()
