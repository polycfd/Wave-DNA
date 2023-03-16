import numpy as np
import matplotlib.pyplot as plt
import sys


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

def _read(filename):
    fo = open(filename, 'r')
    t = []
    p = []
    for line in fo:
        if '#' not in line:
            splitline = line.split()
            t.append(float(splitline[1])*fa)
            p.append(float(splitline[5])/pa)
    fo.close()
    return t,p


optionsfile = 'run.DNA' # default if no command line option (plot *.py optionsfile) is given
if len(sys.argv) > 1:
    optionsfile = sys.argv[1]
fo = open(optionsfile, "r")
for line in fo:
  if "SoundSpeed" in line:
    splitline = line.split()
    c0 = float(splitline[1])
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

la = c0/fa # wavelength
ah = 2.0*c0**2/rh # acceleration magnitude at rh
He = c0**2/(ah*la) # Helmhotz number

filename = "results/horizon.dat"
t,p = _read(filename)
plt.plot(t,p, linewidth=2.0, color='black', linestyle='solid', label=r'$\mathrm{simulation\:linear\:wave}$')


plt.xlabel('$f_{\mathrm{a}}t$', fontsize=fontsize)
plt.ylabel('$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)
plt.legend(loc='upper center', frameon=False, ncol=1, fontsize=fontsize, bbox_to_anchor=[0.5, 2.0])
plt.ylim(0.7,1.1)
plt.xlim(7,15)
plt.xticks([7,8,9,10,11,12,13,14,15], ['$7$','$8$','$9$','$10$','$11$','$12$','$13$','$14$','$15$'], fontsize=fontsize)
plt.yticks([0.7,0.8,0.9,1.0,1.1], ['$0.7$', '$0.8$', '$0.9$', '$1.0$', '$1.1$'], fontsize=fontsize)

plt.grid()
#plt.savefig("../../figures/horizon_smallABH.pdf", bbox_inches = "tight", dpi=400)
#plt.close()
plt.show()

