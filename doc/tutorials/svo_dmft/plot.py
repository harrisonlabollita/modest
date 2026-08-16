import numpy as np

import matplotlib.pyplot as plt

from h5 import HDFArchive

import triqs_modest as tm

with HDFArchive("dmft_results.h5", 'r') as ar:
    mu = ar['mu']
    Sigma_iw = ar['Sigma_iw']

om = np.fromiter(Sigma_iw.mesh, complex).imag
fig, ax = plt.subplots(1, 2, figsize=(6.75, 2.75))
ax[0].plot(om, Sigma_iw[0][0,0].data.imag, 'o-', lw=0.25, ms=4, mec='k')
ax[1].plot(om, Sigma_iw[0][0,0].data.real, 'o-', lw=0.25, ms=4, mec='k')
for a in ax: a.set_xlim(0, 10); a.set_xlabel(r'$\omega_{n}$ (eV)', fontsize=12)
ax[0].set_ylabel(r'$\mathrm{Im}\,\Sigma(i\omega_{n})$ (eV)', fontsize=12)
ax[1].set_ylabel(r'$\mathrm{Re}\,\Sigma(i\omega_{n})$ (eV)', fontsize=12)
fig.tight_layout()
fig.savefig('media/sigma_iw.png', dpi=200)

plt.show()
