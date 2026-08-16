import os

import numpy as np

import matplotlib.pyplot as plt

from h5 import HDFArchive

from triqs.gfs import BlockGf, MeshReFreq
from triqs.lattice.utils import k_space_path

import triqs_modest as tm

DATA = os.path.dirname(__file__)
PLOT = os.path.join(os.path.dirname(__file__), 'media')

TB_FILE = os.path.join(DATA, 'sro')
SW_FILE = os.path.join(DATA, 'Sigma_w.h5')

class BZ:
    A = np.array([[-1.9379489119,  1.9379489110, 6.4363736866999997],
                  [1.93794891190, -1.9379489119,  6.4363736866999997],
                  [1.93794891190,  1.9379489119, -6.4363736866999997]])
    units = 2*np.pi*np.linalg.inv(A).T

bz = BZ()

def Hsoc(lambda_soc):
    lam_xy, lam_z = lambda_soc, lambda_soc
    lam_loc = np.zeros((6,6),dtype=complex)
    lam_loc[0,4] =  1j*lam_xy/2.0
    lam_loc[0,5] =     lam_xy/2.0
    lam_loc[1,2] =  1j*lam_z/2.0
    lam_loc[1,3] = -1j*lam_xy/2.0
    lam_loc[2,3] =    -lam_xy/2.0
    lam_loc[4,5] = -1j*lam_z/2.0
    lam_loc = lam_loc + np.transpose(np.conjugate(lam_loc))
    return lam_loc

eF =  11.3934

G = [0.,0.,0.]; M = [0.5, 0.5, 0.]; X = [0.5, 0.0, 0.]
klabels  = [r'$\Gamma$', r'$\mathrm{X}$', r'$\mathrm{M}$', r'$\Gamma$']
segments = [(G, X), (X, M), (M, G)]
k_vecs, k_plot, k_ticks = k_space_path(segments, num=100, bz=bz)

obe = tm.one_body_elements_from_wannier90(TB_FILE, "NonPolarized", [tm.AtomicOrbs(dim=3)])
obe = tm.add_local_term(tm.extend_to_spin(obe), Hsoc(0.2))

with HDFArchive(SW_FILE) as ar:
    Sigma_dmft = ar['Sw']

om = np.fromiter(Sigma_dmft.mesh,float)
Sigma_phen = Sigma_dmft.copy()
for bl, sig in Sigma_phen:
    Sigma_phen[bl] -= sig(0.0).real # -ReΣ(0)

E = tm.make_embedding(['ud'], [[1,1,1],[1,1,1]], [0,0])
Sigma_C = E.embed([Sigma_phen])

broadening = 1e-7

# A(k,w) along the high-symmetry path
Akw = tm.spectral_function_on_high_symmetry_path(obe, k_vecs, eF, Sigma_C, broadening=broadening)

# A(k,0) on a 2D (kx,ky) grid at kz=0 -> Fermi surface
nk_fs = 101
kx = np.linspace(0, 0.5, nk_fs)
ky = np.linspace(0, 0.5, nk_fs)
KX, KY = np.meshgrid(kx, ky)
k_grid = np.column_stack([KX.ravel(), KY.ravel(), np.zeros(KX.size)])

Akw_fs = tm.spectral_function_on_high_symmetry_path(obe, k_grid, eF, Sigma_C, broadening=0.01)
iw0 = np.argmin(np.abs(om))
A0 = Akw_fs.total[0][:, iw0].reshape(nk_fs, nk_fs)

fig, ax = plt.subplots(1, 3, figsize=(9, 3))

# panel 0: ReΣ(w) - ReΣ(0), xy and xz/yz orbitals only
ax[0].plot(om, Sigma_phen['up_0'].data[:, 0, 0].real, label=r'$xz/yz$', lw=1.2, color='C0')
ax[0].plot(om, Sigma_phen['up_2'].data[:, 0, 0].real, label=r'$xy$', lw=1.2, color='C1')
ax[0].axhline(0, color='k', lw=0.5, ls='dotted')
ax[0].axvline(0, color='k', lw=0.5, ls='dotted')
ax[0].set_xlim(-1, 1); ax[0].set_ylim(-0.5, 0.5)
ax[0].tick_params(which='major', length=4); ax[0].tick_params(which='minor', length=2)
ax[0].set_xlabel(r'$\omega$ $\mathrm{(eV)}$')
ax[0].set_ylabel(r'$\mathrm{Re}\Sigma(\omega) - \mathrm{Re}\Sigma(0)$ $\mathrm{(eV)}$')
ax[0].legend(fontsize=6)
ax[0].set_box_aspect(1)

# panel 1: A(k,w) along high-symmetry path
K, W = np.meshgrid(k_plot, om)
d = ax[1].pcolormesh(K, W, Akw.total[0].T, vmin=0, vmax=8, cmap='magma', shading='auto')
ax[1].axhline(0, color='xkcd:yellow', lw=0.5, ls='dotted')
ax[1].set_ylim(-0.6, +0.5)
ax[1].set_ylabel(r'$\omega$ $\mathrm{(eV)}$')
ax[1].set_xticks(k_ticks, labels=klabels)
ax[1].tick_params(axis='x', which='both', length=0)
ax[1].tick_params(axis='y', direction='out', length=3)
ax[1].tick_params(axis='y', which='major', direction='out', length=4)
ax[1].tick_params(axis='y', which='minor', direction='out', length=0)
cax1 = ax[1].inset_axes([0.0, 1.05, 1.0, 0.05])
cbar = fig.colorbar(d, cax=cax1, orientation='horizontal')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.tick_params(which='both', length=0)
cbar.ax.tick_params(which='major', length=2, direction='out')
ax[1].set_box_aspect(1)

# panel 2: A(k,0) on the (kx,ky) grid -- Fermi surface
d = ax[2].pcolormesh(KX, KY, A0, cmap='magma', shading='auto', vmin=0)
ax[2].set_xlabel(r'$k_x$')
ax[2].set_ylabel(r'$k_y$')
ax[2].set_box_aspect(1)
cax2 = ax[2].inset_axes([0.0, 1.05, 1.0, 0.05])
cbar = fig.colorbar(d, cax=cax2, orientation='horizontal')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.tick_params(which='both', length=0)
cbar.ax.tick_params(which='major', length=2, direction='out')
ax[2].set_xticks([0, 0.5], labels=[r'$0$', r'$\pi$']); ax[2].set_yticks([0., 0.5], labels=[r'$0$', r'$\pi$'])

for a, label in zip(ax, ['(a)', '(b)', '(c)']):
    a.text(-0.25, 1.15, label, transform=a.transAxes, ha='left', va='top',
           fontsize=10)

plt.subplots_adjust(wspace=0.45)
out = os.path.join(PLOT, 'sro_sigma_akw_fs.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
