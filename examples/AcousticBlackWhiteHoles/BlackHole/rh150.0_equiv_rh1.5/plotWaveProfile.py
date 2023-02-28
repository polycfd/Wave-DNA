import numpy as np
import matplotlib.pyplot as plt

#===================================================================
preamble='\\usepackage{amsmath}\n\\usepackage{amsfonts}'
plt.rcParams['mathtext.fontset']='cm' 
fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6,2))
useMathText=True
fontsize = 20
ax.xaxis.label.set_size(fontsize)
ax.yaxis.label.set_size(fontsize)
plt.setp(ax.get_xticklabels(), fontsize=fontsize)
plt.setp(ax.get_yticklabels(), fontsize=fontsize)
#===================================================================

rh = 1.0 # default reference length scale if case is not ABH or AWH

fo = open("run.DNA", "r")
for line in fo:
  if "ExcitationFrequency" in line:
    splitline = line.split()
    fa = float(splitline[1])
  if "PressureAmplitude" in line:
    splitline = line.split()
    pa = float(splitline[1])
  if "HorizonRadius" in line:
    splitline = line.split()
    rh = float(splitline[1])
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
                    x.append(float(splitline[1])/rh)
                    p.append(float(splitline[2])/pa)
                    line = next(fo)
    fo.close()
    return x,p

filename = "results/fields.dat"

x,p = _read(filename, 40000)
plt.plot(x,p, linewidth=3.0, color='red', label='$f_{\mathrm{a}}t=10$')
x,p = _read(filename, 60000)
plt.plot(x,p, linewidth=1.5, color='black', label='$f_{\mathrm{a}}t=15$')

plt.xlabel('$r/r_{\mathrm{h}}$', fontsize=fontsize)
plt.ylabel('$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)
plt.ylim(-1.5,1.5)
plt.xlim(0.9,1.1)
plt.xticks([0.90, 0.95, 1.00, 1.05, 1.10], ['$0.90$', '$0.95$', '$1.00$', '$1.05$', '$1.10$'], fontsize=fontsize)
plt.yticks([-1, 0, 1], ['$-1$', '$0$', '$1$'], fontsize=fontsize)

plt.legend(loc='upper center', frameon=False, ncol=2, fontsize=fontsize, bbox_to_anchor=[0.5, 1.4])
plt.grid()
#plt.savefig("../../figures/pr_equivalentABH.pdf", bbox_inches = "tight", dpi=400)
#plt.close()
plt.show()
