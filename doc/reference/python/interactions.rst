.. _reference/python/interactions:

Interaction Hamiltonians
************************

ModEST provides factory functions that build a local interaction
Hamiltonian as a TRIQS many-body operator using the "flavors" of the
impurity problem defined by the embedding. The generic form is

.. math::

   H_{\mathrm{int}} =
   \frac{1}{2} \sum_{\substack{\alpha\beta\gamma\delta \\ \sigma\sigma'}}
   U_{\alpha\beta\gamma\delta}\;
   c^{\dagger}_{\alpha\sigma}\, c^{\dagger}_{\beta\sigma'}\,
   c_{\delta\sigma'}\, c_{\gamma\sigma},

where :math:`\alpha,\beta,\gamma,\delta` index the orbitals (or flavors) of the
impurity, :math:`\sigma,\sigma'` are spin (or block-diagonal) indices,
and :math:`U_{\alpha\beta\gamma\delta}` is the matrix of two-particle
interactions. The factories below differ in how
:math:`U_{\alpha\beta\gamma\delta}` is parametrized — density–density
only, Kanamori, or full rotationally-invariant Slater — and return a
``triqs.operators.Operator`` ready to be handed to a TRIQS impurity
solver.

.. autosummary::

   triqs_modest.hamiltonians.make_density_density
   triqs_modest.hamiltonians.make_kanamori
   triqs_modest.hamiltonians.make_slater

Coulomb tensors
---------------

``make_slater`` takes the four-index Coulomb tensor as an argument, which 
can be parametrized by :math:`U` and :math:`J`. A typical use is::

   U_shell = tm.slater_tensor(obe.C_space, atom=0, U_int=U, J_hund=J)   # (2l+1)^4
   U_local = tm.to_local_basis(U_shell, obe.C_space, atom=0)            # dim^4
   h_int   = tm.make_slater(E, 0, U_local)

The first line constructs the tensor for the full :math:`2l+1` shell of the first
atom in the C space, in the orbital basis of the DFT code: it is built in the
spherical harmonic basis and rotated by the :math:`Y_l^m \rightarrow` DFT rotation,
both read from the local space. The second line restricts it to the orbitals the
correlated space spans and rotates it into the local basis. The third constructs
the interaction Hamiltonian.

If the correlated space covers only part of the shell (a t2g or :math:`e_g`
subset), pass the shell positions of the orbitals you want sliced out to
``to_local_basis``, in the orbital ordering of the DFT code — for VASP's
:math:`l=2` ordering (xy, yz, z², xz, x²-y²) that is ``orbs=[0,1,3]`` for t2g and
``orbs=[2,4]`` for :math:`e_g`.

Another scenario is when the impurity only covers part of the correlated space —
for example an impurity made with ``split_imp``. Nothing in the three lines above
changes: ``make_slater`` takes the embedding ``E``, so it can cut the impurity's
block out of the :math:`\dim(\mathcal{C})^4` tensor, provided its orbitals are
contiguous in :math:`\mathcal{C}`::

   E     = tm.make_embedding(obe.C_space).split_imp(0, [0, 1, 2])
   h_int = tm.make_slater(E, 0, U_local)   # imp 0's 3x3x3x3 block is cut out

If they are not contiguous, restrict the tensor yourself and pass the
impurity-sized one.

.. autosummary::

   triqs_modest.hamiltonians.slater_tensor
   triqs_modest.hamiltonians.to_local_basis
   triqs_modest.hamiltonians.U_matrix_slater_spherical
   triqs_modest.hamiltonians.U_matrix_slater_local
   triqs_modest.hamiltonians.U_matrix_kanamori
   triqs_modest.hamiltonians.rotate_U_matrix_slater
