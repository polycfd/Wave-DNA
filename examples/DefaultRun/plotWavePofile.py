import numpy as np
import matplotlib.pyplot as plt


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

# Default settings
pa = 1.0e3
fa = 1.0e5
c0 = 1500.0

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

x,p = _read(filename, 10000)
plt.plot(x,p, linewidth=2.0, color='black')

plt.xlabel('$x/\lambda_0$', fontsize=fontsize)
plt.ylabel('$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize)
#plt.legend(loc='upper right', frameon=True, ncol=1, fontsize=fontsize)

plt.xlim(0,10)
plt.ylim(-1.1,1.1)

plt.grid()
plt.show()
#plt.savefig("../figures/DefaultRun.pdf", bbox_inches = "tight")
#plt.close()

