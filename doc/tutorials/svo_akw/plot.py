import numpy as np

import matplotlib.pyplot as plt

from h5 import HDFArchive

with HDFArchive('svo_akw.h5', 'r') as ar:
    Akw = ar['Akw']

n_k, n_w = Akw.shape
w_min, w_max = -1, +1
om = np.linspace(w_min, w_max, n_w)

k_info = [(1, r'$\Gamma$'), (114, 'X'), (227, 'M'), (387, r'$\Gamma$'), (500, 'Z')]
k_ticks, k_labels = map(list, zip(*k_info))
k_plot = np.linspace(min(k_ticks), max(k_ticks), n_k)

fig, ax = plt.subplots(figsize=(4.5, 3.5))
K, W = np.meshgrid(k_plot, om)
d = ax.pcolormesh(K, W, Akw.T, vmin=0, vmax=10, cmap='terrain', shading='auto')
ax.axhline(0, color='w', lw=0.5, ls='dotted')
ax.set_ylabel(r'$\omega$ $\mathrm{(eV)}$')
ax.set_xticks(k_ticks, labels=k_labels)
ax.tick_params(axis='x', which='both', length=0)
ax.tick_params(axis='y', which='major', direction='out', length=4)
ax.tick_params(axis='y', which='minor', direction='out', length=0)

cax = ax.inset_axes([1.03, 0.0, 0.04, 1.0])
cbar = fig.colorbar(d, cax=cax, orientation='vertical')
cbar.ax.tick_params(which='major', length=2, direction='out')
cbar.set_label(r'$A(\mathbf{k},\omega)$')

fig.tight_layout()
fig.savefig('media/svo_akw.png', dpi=200, bbox_inches='tight')

plt.show()
