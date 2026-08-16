import numpy as np
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.gfs import BlockGf, MeshReFreq
import triqs_modest as tm

mu, Sigma_w = None, None
if mpi.is_master_node():
    with HDFArchive('dmft_results.h5') as ar:
        mu = ar['mu']
        Sigma_w = ar['Sigma_w_m_dc']

mu, Sigma_w = mpi.bcast(mu), mpi.bcast(Sigma_w)

_, obe = tm.one_body_elements_from_dft_converter("svo_dft_data.h5")
obe = tm.one_body_elements_on_high_symmetry_path("svo_dft_data.h5", obe)

E = tm.make_embedding(obe.C_space)

Akw = tm.spectral_function_on_high_symmetry_path(obe, mu, E.embed([Sigma_w]), broadening=0.005)

if mpi.is_master_node():
    with HDFArchive('svo_akw.h5') as ar:
        ar['Akw'] = Akw.total[0]
