import numpy as np

import matplotlib.pyplot as plt

from h5 import HDFArchive

from triqs.gfs import BlockGf, MeshImFreq, MeshReFreq

import triqs_modest as tm

with HDFArchive("dmft_results.h5", 'r') as ar:
    mu = ar['mu']
    Sigma_iw = ar['Sigma_iw']

U, J      = 4.5, 0.68
nominal_d = 1.0
sig_dc    = (U - 2*J)*(nominal_d - 0.5)

w_min, w_max = -1, +1
mesh = MeshReFreq(n_w=10001, w_min=w_min, w_max=w_max)
gf_struct = [(bl, g.target_shape[0]) for (bl, g) in Sigma_iw]
Sigma_w = BlockGf(mesh=mesh, gf_struct=gf_struct)
for bl, gf in Sigma_w: gf.set_from_pade(Sigma_iw[bl], n_points=301, freq_offset=0.0)
Sigma_w_m_dc = Sigma_w.copy()
for bl, gf in Sigma_w_m_dc: gf -= sig_dc


om = np.fromiter(Sigma_w.mesh, float)
fig,ax = plt.subplots(1, 2, figsize=(6.75, 2.75))
ax[0].plot(om, -Sigma_w[0][0,0].data.imag, '-', lw=1.25)
ax[1].plot(om, Sigma_w[0][0,0].data.real, '-', lw=1.25)
for a in ax: a.set_xlabel(r'$\omega$ (eV)', fontsize=12)
ax[0].set_ylabel(r'$\mathrm{Im}\Sigma(\omega)$ (eV)', fontsize=12)
ax[1].set_ylabel(r'$\mathrm{Re}\Sigma(\omega)$ (eV)', fontsize=12)
fig.tight_layout()
fig.savefig('media/sigma_w.png', dpi=200)

with HDFArchive('dmft_results.h5', 'a') as ar:
    ar['Sigma_w']  = Sigma_w
    ar['Sigma_w_m_dc']  = Sigma_w_m_dc

get_conductivity = lambda Phi: tm.optical_conductivity(Phi, beta=40.0, contribution='inter').sigma['xx']
with HDFArchive('svo_dmft_optics.h5', 'r') as ar:
    Om = ar['Omega']
    sig_dmft = get_conductivity(ar['dmft']['Phi'])
    sig_dft  = get_conductivity(ar['dft']['Phi'])
plt.figure(figsize=(3.375,3))
plt.plot(Om, sig_dft, 'k-', label='DFT')
plt.plot(Om, sig_dmft, '-', color='orange', label='DMFT')
plt.axvline(0.07, ls='--', lw=1, color='0.5', label='exp. (Ahn et al.)')  # interband transition ~70 meV
plt.xlabel(r'$\Omega$ (eV)'); plt.ylabel(r'$\sigma_{1}^{\mathrm{inter}}$')
plt.xlim(0, 0.3); plt.ylim(0, 3.1)
plt.legend(loc='best', fontsize=10, frameon=False)
plt.tight_layout()
plt.savefig('media/sigma_optics.png', dpi=200)

plt.show()
