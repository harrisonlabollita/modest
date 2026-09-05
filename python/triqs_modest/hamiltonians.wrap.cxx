
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// ==================== module functions ====================

// U_matrix_kanamori
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](long n_orb, double U_int, double U_prime, double J_hund) { return triqs::U_matrix_kanamori(n_orb, U_int, U_prime, J_hund); },
              "n_orb", "U_int", "U_prime", "J_hund")};

// U_matrix_slater_local
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](long l, nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> s2l,
      double U_int, double J_hund) { return triqs::U_matrix_slater_local(l, s2l, U_int, J_hund); },
   "l", "s2l", "U_int", "J_hund")};

// U_matrix_slater_spherical
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](long l, double U_int, double J_hund) { return triqs::U_matrix_slater_spherical(l, U_int, J_hund); }, "l", "U_int", "J_hund")};

// make_density_density
static auto const _c2py_fun_3 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](const triqs::modest::embedding &E, long imp_idx, double U_int, double U_prime,
                                         double J_hund) { return triqs::make_density_density(E, imp_idx, U_int, U_prime, J_hund); },
                                      "E", "imp_idx", "U_int", "U_prime", "J_hund")};

// make_kanamori
static auto const _c2py_fun_4 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const triqs::modest::embedding &E, long imp_idx, double U_int, double U_prime, double J_hund, bool spin_flip,
                 bool pair_hopping) { return triqs::make_kanamori(E, imp_idx, U_int, U_prime, J_hund, spin_flip, pair_hopping); },
              "E", "imp_idx", "U_int", "U_prime", "J_hund", "spin_flip"_a = true, "pair_hopping"_a = true)};

// make_slater
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const triqs::modest::embedding &E, long imp_idx,
      const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>
         &U_tensor) { return triqs::make_slater(E, imp_idx, U_tensor); },
   "E", "imp_idx", "U_tensor")};

// rotate_U_matrix_slater
static auto const _c2py_fun_6 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>
         &U_tensor,
      nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> M) {
     return triqs::rotate_U_matrix_slater(U_tensor, M);
   },
   "U_tensor", "M")};

// slater_tensor
static auto const _c2py_fun_7 = c2py::dispatcher_f_kw_t{c2py::cfun([](const triqs::modest::local_space &C_space, long atom, double U_int,
                                                                      double J_hund) { return triqs::slater_tensor(C_space, atom, U_int, J_hund); },
                                                                   "C_space", "atom", "U_int", "J_hund")};

// to_local_basis
static auto const _c2py_fun_8 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>
         &U_tensor,
      const triqs::modest::local_space &C_space, long atom,
      const std::optional<std::vector<long>> &orbs) { return triqs::to_local_basis(U_tensor, C_space, atom, orbs); },
   "U_tensor", "C_space", "atom", "orbs"_a = std::optional<std::vector<long>>{})};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(
   R"DOC(
Construct the U and U' interaction matrices for the Hubbard-Kanamori model.

Parameters
----------
n_orb : {par_0}
   Number of orbitals.
U_int : {par_1}
   Screened Hubbard interaction.
U_prime : {par_2}
   :math:`U'` (typically :math:`U' = U - 2J`).
J_hund : {par_3}
   Hund's coupling.

Returns
-------
{ret_0}
   pair of interaction matrices (U, U').
)DOC",
   {{c2py::python_typename<long>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}},
   {c2py::python_typename<
      std::pair<nda::basic_array<double, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>,
                nda::basic_array<double, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(
   R"DOC(
Construct a four-index Coulomb tensor in a specific orbital basis.

Parameters
----------
l : {par_0}
   Angular quantum number.
s2l : {par_1}
   Spherical to local basis rotation.
U_int : {par_2}
   Screened Hubbard interaction.
J_hund : {par_3}
   Hund's coupling.

Returns
-------
{ret_0}
   Coulomb tensor.
)DOC",
   {{c2py::python_typename<long>()},
    {c2py::python_typename<
       nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()},
    {c2py::python_typename<double>()},
    {c2py::python_typename<double>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_2 = _c2py_fun_2.doc(
   R"DOC(
Construct a four-index Coulomb tensor in the basis of spherical harmonics.

We typically construct the four-index Coulomb tensor in the basis of spherical harmonics,

.. math::

   U_{m_{1}m_{2}m_{3}m_{4}}^{\mathrm{spherical}} = \sum_{k=0}^{2l} F_{k} \alpha (l, k, m_{1}, m_{2}, m_{3}, m_{4}),

where :math:`F_{k}` are radial Slater integrals and :math:`\alpha(l, k, m_{1}, m_{2}, m_{3}, m_{4})` denote angular
Racah-Wigner numbers for a spherically symmetric interaction tensor.

Parameters
----------
l : {par_0}
   Angular quantum number.
U_int : {par_1}
   Screened Hubbard interaction.
J_hund : {par_2}
   Hund's coupling.

Returns
-------
{ret_0}
   Coulomb tensor.
)DOC",
   {{c2py::python_typename<long>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}},
   {c2py::python_typename<nda::basic_array<double, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Construct a density-density interaction Hamiltonian for one impurity.

Create a density-density Hamiltonian

.. math::

   H_{\mathrm{int}} = \frac{1}{2} \sum_{(i\sigma)\neq(j\sigma^{\prime})} U_{ij}^{\sigma\sigma^{\prime}}n_{i\sigma}
   n_{j\sigma^{\prime}}.

Equivalent to ``triqs::make_kanamori`` with the spin-flip and pair-hopping terms switched off.

Parameters
----------
E : {par_0}
   The embedding.
imp_idx : {par_1}
   Index of the impurity.
U_int : {par_2}
   Hubbard :math:`U`.
U_prime : {par_3}
   :math:`U'` (typically :math:`U' = U - 2J`).
J_hund : {par_4}
   Kanamori :math:`J`.

Returns
-------
{ret_0}
   Many-body operator representing the Hamiltonian.
)DOC",
                                                {{c2py::python_typename<const triqs::modest::embedding &>()},
                                                 {c2py::python_typename<long>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<double>()}},
                                                {c2py::python_typename<triqs::operators::many_body_operator>()});
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC(
Construct a Hubbard-Kanamori Hamiltonian for one impurity.

Create a Hubbard-Kanamori Hamiltonian using the density-density, spin-flip, and pair-hopping interactions,

.. math::

   H_{\mathrm{int}} = \frac{1}{2} \sum_{(i\sigma)\neq(j\sigma^{\prime})} U_{ij}^{\sigma\sigma^{\prime}}n_{i\sigma}
   n_{j\sigma^{\prime}} - \sum_{i\neq j}Jc_{i\uparrow}^{\dagger}c_{i\downarrow}c_{j\downarrow}^{\dagger}
   c_{j\uparrow} + \sum_{i\neq j} J c_{i\uparrow}^{\dagger}c_{i\downarrow}^{\dagger}c_{j\downarrow}c_{j\uparrow}.

The operator names follow the block structure of impurity `imp_idx`, i.e. `(tau_gamma, orbital)`.

Parameters
----------
E : {par_0}
   The embedding.
imp_idx : {par_1}
   Index of the impurity.
U_int : {par_2}
   Hubbard :math:`U`.
U_prime : {par_3}
   :math:`U'` (typically :math:`U' = U - 2J`).
J_hund : {par_4}
   Kanamori :math:`J`.
spin_flip : {par_5}
   Spin flip term.
pair_hopping : {par_6}
   Pair-hopping term.

Returns
-------
{ret_0}
   Many-body operator representing the Hamiltonian.
)DOC",
                                                {{c2py::python_typename<const triqs::modest::embedding &>()},
                                                 {c2py::python_typename<long>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<double>()},
                                                 {c2py::python_typename<bool>()},
                                                 {c2py::python_typename<bool>()}},
                                                {c2py::python_typename<triqs::operators::many_body_operator>()});
static const auto _c2py_doc_5 =
   _c2py_fun_5.doc(R"DOC(
Construct a Slater Hamiltonian for one impurity from a Coulomb tensor.

Create a Slater Hamiltonian using fully rotationally-invariant four-index interactions:

.. math::

   H_{\mathrm{int}} = \frac{1}{2} \sum_{ijkl, \sigma\sigma^{\prime}} U_{ijkl}c^{\dagger}_{i\sigma}
   c^{\dagger}_{j\sigma^{\prime}}c_{l\sigma^{\prime}}c_{k\sigma}.

`U_tensor` may be given at either of two sizes:

* sized for impurity `imp_idx` itself, the usual case when the impurity is a whole atomic shell;
* sized for the whole correlated space :math:`\mathcal{C}`, in which case this impurity's block is cut out
  through the embedding. This is what ``triqs::to_local_basis`` returns when the impurity is a part of
  :math:`\mathcal{C}` split off with `split_imp`, e.g. a t2g impurity from a full d shell. It requires the
  impurity's orbitals to form one contiguous block of :math:`\mathcal{C}`; if they do not, restrict the
  tensor yourself and pass the impurity-sized one.

Parameters
----------
E : {par_0}
   The embedding.
imp_idx : {par_1}
   Index of the impurity.
U_tensor : {par_2}
   Coulomb tensor in the local basis, of shape :math:`n_{\mathrm{orb}}^4` or
   :math:`\mathrm{dim}(\mathcal{C})^4`, e.g. built with ``triqs::slater_tensor`` followed by
   ``triqs::to_local_basis.``

Returns
-------
{ret_0}
   Many-body operator representing the Hamiltonian.
)DOC",
                   {{c2py::python_typename<const triqs::modest::embedding &>()},
                    {c2py::python_typename<long>()},
                    {c2py::python_typename<const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A',
                                                                  nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> &>()}},
                   {c2py::python_typename<triqs::operators::many_body_operator>()});
static const auto _c2py_doc_6 = _c2py_fun_6.doc(
   R"DOC(
Rotate a four-index Coulomb tensor into a new orbital basis.

`M` transforms the orbitals, :math:`|i\rangle_{\mathrm{new}} = \sum_{j} M_{ij}|j\rangle_{\mathrm{old}}`,
so the operators pick up the conjugate, :math:`c^{\mathrm{new}}_{i} = \sum_{j} M^{*}_{ij}c^{\mathrm{old}}_{j}`,
and the Coulomb tensor transforms as

.. math::

   U'_{iknp} = \sum_{jqmo} M^{*}_{ij} M^{*}_{kq} U_{jqmo} M_{nm} M_{po}.

Parameters
----------
U_tensor : {par_0}
   Coulomb tensor in the old basis.
M : {par_1}
   Orbital change-of-basis matrix (new :math:`\leftarrow` old).

Returns
-------
{ret_0}
   Coulomb tensor in the new basis.
)DOC",
   {{c2py::python_typename<
       const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> &>()},
    {c2py::python_typename<
       nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_7 = _c2py_fun_7.doc(
   R"DOC(
Build the Slater Coulomb tensor for one atomic shell, in the basis of its DFT code.

The angular momentum :math:`l` and the :math:`Y_l^m \rightarrow` DFT rotation are both read from
`C_space`, so neither has to be supplied (and neither can be inferred from the number of correlated orbitals:
a t2g or :math:`e_g` shell has `dim` :math:`< 2l+1`).

The returned tensor always spans the **full** :math:`(2l+1)` shell. If the correlated space covers only part of
it, restrict the tensor afterwards with ``triqs::to_local_basis.``

Parameters
----------
C_space : {par_0}
   The local :math:`\mathcal{C}` space, e.g. `obe.C_space`.
atom : {par_1}
   Index of the atomic shell within `C_space`.
U_int : {par_2}
   Screened Hubbard interaction.
J_hund : {par_3}
   Hund's coupling.

Returns
-------
{ret_0}
   Coulomb tensor of shape :math:`(2l+1)^4` in the DFT orbital basis.
)DOC",
   {{c2py::python_typename<const triqs::modest::local_space &>()},
    {c2py::python_typename<long>()},
    {c2py::python_typename<double>()},
    {c2py::python_typename<double>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_8 = _c2py_fun_8.doc(
   R"DOC(
Restrict a shell Coulomb tensor to the correlated orbitals and rotate it into the local basis.

Performs the two steps that take a full-shell tensor to one that can be used with the impurity
Green's functions:

1. Select the `orbs` orbitals of the shell that the correlated space actually spans, giving a
   :math:`\mathrm{dim}^4` tensor.
2. Change basis to match the one-body elements, which are rotated as :math:`P \leftarrow R^{\dagger}P` with
   :math:`R` the stored DFT-to-local rotation of the shell. The operators therefore transform with
   :math:`R^{\dagger}`, so the *orbital* matrix passed to ``triqs::rotate_U_matrix_slater`` is its
   conjugate, :math:`R^{T}`.

`orbs` indexes positions within the :math:`(2l+1)` shell, in the orbital ordering of the DFT code, and is
ordered: entry `i` names the shell orbital that is correlated orbital `i`. The selection need not be
contiguous -- for VASP's :math:`l=2` ordering (xy, yz, z², xz, x²-y²), t2g is `{0,1,3}` and :math:`e_g` is
`{2,4}`.

If the correlated space spans the whole shell (`dim == 2l+1`), `orbs` may be omitted. Otherwise it is
required: which orbitals the projectors span is not recorded by the DFT converters and cannot be inferred.

:math:`R` is read from `C_space`, i.e. it is the rotation the loader applied. `rotate_local_basis` does not
record itself there, so if you rotated the one-body elements yourself, rotate the tensor yourself too, with
``triqs::rotate_U_matrix_slater`` and the transpose of the rotation you passed.

Parameters
----------
U_tensor : {par_0}
   Coulomb tensor spanning the full shell, e.g. from ``triqs::slater_tensor.``
C_space : {par_1}
   The local :math:`\mathcal{C}` space.
atom : {par_2}
   Index of the atomic shell within `C_space`.
orbs : {par_3}
   Shell positions of the correlated orbitals; defaults to the whole shell.

Returns
-------
{ret_0}
   Coulomb tensor of shape :math:`\mathrm{dim}^4` in the local basis.
)DOC",
   {{c2py::python_typename<
       const nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> &>()},
    {c2py::python_typename<const triqs::modest::local_space &>()},
    {c2py::python_typename<long>()},
    {c2py::python_typename<const std::optional<std::vector<long>> &>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 4, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"U_matrix_kanamori", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"U_matrix_slater_local", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"U_matrix_slater_spherical", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"make_density_density", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"make_kanamori", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"make_slater", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"rotate_U_matrix_slater", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"slater_tensor", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"to_local_basis", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "hamiltonians",                                                                                                  /* name of module */
   R"RAWDOC(Functions and utilities to construct interaction Hamiltonians using TRIQS many-body operators.)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_hamiltonians() {

  if (not c2py::check_python_version("hamiltonians")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)

#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
