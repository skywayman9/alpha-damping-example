from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')

plt.rcParams.update({'font.size': 12})


x,y,d,p,u,v=np.loadtxt('Rslt0047.plt', delimiter=None, unpack=True,skiprows=3)
# x,y,d,p,u,v=np.loadtxt('soln.txt', delimiter=None, unpack=True)


g=320
k=640


x = x.reshape(g,k)
y = y.reshape(g,k)
d = d.reshape(g,k)
p = p.reshape(g,k)
u = u.reshape(g,k)
v = v.reshape(g,k)

print(np.max(d))
plt.xlim(0.4,1.0)
plt.ylim(0,0.5)
plt.contour(x,y,d,43,linewidths=0.5,colors=('k'))
plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()

fig1.savefig('FD_500.pdf', dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()