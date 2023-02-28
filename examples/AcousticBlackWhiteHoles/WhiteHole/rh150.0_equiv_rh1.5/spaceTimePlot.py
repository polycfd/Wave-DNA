import numpy as np
import matplotlib.pyplot as plt

#===================================================================
preamble='\\usepackage{amsmath}\n\\usepackage{amsfonts}'
plt.rcParams['mathtext.fontset']='cm' 
fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6,3))
useMathText=True
fontsize = 20
ax.xaxis.label.set_size(fontsize)
ax.yaxis.label.set_size(fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
#===================================================================


rh = 1.0 # default reference length scale if case is not ABH or AWH

fo = open("run.DNA", "r")
for line in fo:
  if "TimeStepSize" in line:
    splitline = line.split()
    dt = float(splitline[1])
  if "NPoints" in line:
    splitline = line.split()
    Npoints = int(splitline[1])
  if "BoundaryMotionStartTime" in line:
    splitline = line.split()
    StartTime = 0.0 #float(splitline[1])
  if "FixedBoundaryPosition" in line:
    splitline = line.split()
    Xmax = float(splitline[1])
  if "WriteFrequency" in line:
    splitline = line.split()
    WriteFrequency = float(splitline[1])
  if "TimeEnd" in line:
    splitline = line.split()
    EndTime = float(splitline[1])
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





TimeRange = EndTime - StartTime
Nsets = int(TimeRange/(dt*WriteFrequency)-1)

t = np.linspace(0.0, TimeRange*fa, num=Nsets)
x = np.linspace(0.0, Xmax, num=Npoints)
X, T = np.meshgrid(x, t)
np.transpose(T)
P = np.zeros((Nsets, Npoints))

fo = open('results/fields.dat', "r")
Set = 0
for line in fo:
  if "timeStep" in line and Set<Nsets:
    line = next(fo)
    splitline = line.split()
    R = float(splitline[1])
    X[Set,:] = np.linspace(R, Xmax, num=Npoints)/rh
    line = next(fo)
    line = next(fo)
    i = 0
    while "EOS" not in line:
      splitline = line.split()
      P[Set, i] = float(splitline[2])/pa
      line = next(fo)
      i = i+1
    Set = Set+1
fo.close()

cmap = plt.get_cmap('Blues_r')
im = ax.pcolormesh(X[:,0:-1:1], T[:,0:-1:1], P[:,0:-1:1], cmap=cmap, rasterized=True)

cb = fig.colorbar(im, ax=ax, orientation='vertical')
cb.ax.tick_params(labelsize=fontsize)

cb.ax.set_title(r'$p_1/\Delta p_{\mathrm{a}}$', fontsize=fontsize, pad=15)
plt.xlabel(r'$r/r_{\mathrm{h}}$')
plt.ylabel(r'$f_{\mathrm{a}}t$')

plt.xlim((0.9,1.1))
plt.ylim((0, 15.0))
plt.xticks([0.90, 0.95, 1.00, 1.05, 1.10], ['$0.90$', '$0.95$', '$1.00$', '$1.05$', '$1.10$'], fontsize=fontsize)
plt.yticks([0, 5, 10, 15], ['$0$', '$5$', '$10$', '$15$'], fontsize=fontsize)
cb.set_ticks([-1, 0, 1])
cb.set_ticklabels(['$-1$', '$0$', '$1$'])
im.set_clim(-1,1)

plt.legend(fontsize=fontsize, loc='lower right', ncol=1, frameon=False)
ax.set_facecolor('grey')
plt.grid(color='black')
#plt.savefig("../../figures/rt_equivalentAWH.pdf", bbox_inches = "tight", dpi=400)
#plt.close()
plt.show()
  
