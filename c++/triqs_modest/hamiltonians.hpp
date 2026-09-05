// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/gfs.hpp>
#include <triqs/operators.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include "./local_space.hpp"
#include "./embedding.hpp"

using namespace triqs::gfs;

namespace triqs {

  /** @name Coulomb tensor utilities
   * Utility functions for creating and working with the four index Coulomb tensor.
   * @{
   */

  /**
   * @ingroup hamiltonian
   * @brief Construct the U and U' interaction matrices for the Hubbard-Kanamori model.
   *
   * @param n_orb Number of orbitals.
   * @param U_int Screened Hubbard interaction.
   * @param U_prime \f$ U' \f$ (typically \f$ U' = U - 2J \f$).
   * @param J_hund Hund's coupling.
   * @return pair of interaction matrices (U, U').
   */
  std::pair<nda::matrix<double>, nda::matrix<double>> U_matrix_kanamori(long n_orb, double U_int, double U_prime, double J_hund);

  /**
   * @ingroup hamiltonian
   * @brief Construct a four-index Coulomb tensor in the basis of spherical harmonics.
   *
   * @details We typically construct the four-index Coulomb tensor in the basis of spherical harmonics,
   * \f[
   *   U_{m_{1}m_{2}m_{3}m_{4}}^{\mathrm{spherical}} = \sum_{k=0}^{2l} F_{k} \alpha (l, k, m_{1}, m_{2}, m_{3}, m_{4}),
   * \f]
   * where \f$ F_{k} \f$ are radial Slater integrals and \f$ \alpha(l, k, m_{1}, m_{2}, m_{3}, m_{4}) \f$ denote angular
   * Racah-Wigner numbers for a spherically symmetric interaction tensor.
   *
   * @param l Angular quantum number.
   * @param U_int Screened Hubbard interaction.
   * @param J_hund Hund's coupling.
   * @return Coulomb tensor.
   */
  nda::array<double, 4> U_matrix_slater_spherical(long l, double U_int, double J_hund);

  /**
   * @ingroup hamiltonian
   * @brief Construct a four-index Coulomb tensor in a specific orbital basis.
   *
   * @param l Angular quantum number.
   * @param s2l Spherical to local basis rotation.
   * @param U_int Screened Hubbard interaction.
   * @param J_hund Hund's coupling.
   * @return Coulomb tensor.
   */
  nda::array<dcomplex, 4> U_matrix_slater_local(long l, nda::matrix<dcomplex> s2l, double U_int, double J_hund);

  /**
   * @ingroup hamiltonian
   * @brief Rotate a four-index Coulomb tensor into a new orbital basis.
   *
   * @details `M` transforms the orbitals, \f$ |i\rangle_{\mathrm{new}} = \sum_{j} M_{ij}|j\rangle_{\mathrm{old}} \f$,
   * so the operators pick up the conjugate, \f$ c^{\mathrm{new}}_{i} = \sum_{j} M^{*}_{ij}c^{\mathrm{old}}_{j} \f$,
   * and the Coulomb tensor transforms as
   * \f[
   *   U'_{iknp} = \sum_{jqmo} M^{*}_{ij} M^{*}_{kq} U_{jqmo} M_{nm} M_{po}.
   * \f]
   *
   * @param U_tensor Coulomb tensor in the old basis.
   * @param M Orbital change-of-basis matrix (new \f$ \leftarrow \f$ old).
   * @return Coulomb tensor in the new basis.
   */
  nda::array<dcomplex, 4> rotate_U_matrix_slater(nda::array<dcomplex, 4> const &U_tensor, nda::matrix<dcomplex> M);

  /**
   * @ingroup hamiltonian
   * @brief Build the Slater Coulomb tensor for one atomic shell, in the basis of its DFT code.
   *
   * @details The angular momentum \f$ l \f$ and the \f$ Y_l^m \rightarrow \f$ DFT rotation are both read from
   * `C_space`, so neither has to be supplied (and neither can be inferred from the number of correlated orbitals:
   * a t2g or \f$ e_g \f$ shell has `dim` \f$ < 2l+1 \f$).
   *
   * The returned tensor always spans the **full** \f$ (2l+1) \f$ shell. If the correlated space covers only part of
   * it, restrict the tensor afterwards with @ref triqs::to_local_basis.
   *
   * @param C_space The local \f$ \mathcal{C} \f$ space, e.g. `obe.C_space`.
   * @param atom Index of the atomic shell within `C_space`.
   * @param U_int Screened Hubbard interaction.
   * @param J_hund Hund's coupling.
   * @return Coulomb tensor of shape \f$ (2l+1)^4 \f$ in the DFT orbital basis.
   */
  nda::array<dcomplex, 4> slater_tensor(modest::local_space const &C_space, long atom, double U_int, double J_hund);

  /**
   * @ingroup hamiltonian
   * @brief Restrict a shell Coulomb tensor to the correlated orbitals and rotate it into the local basis.
   *
   * @details Performs the two steps that take a full-shell tensor to one that can be used with the impurity
   * Green's functions:
   *
   * 1. Select the `orbs` orbitals of the shell that the correlated space actually spans, giving a
   *    \f$ \mathrm{dim}^4 \f$ tensor.
   * 2. Change basis to match the one-body elements, which are rotated as \f$ P \leftarrow R^{\dagger}P \f$ with
   *    \f$ R \f$ the stored DFT-to-local rotation of the shell. The operators therefore transform with
   *    \f$ R^{\dagger} \f$, so the *orbital* matrix passed to @ref triqs::rotate_U_matrix_slater is its
   *    conjugate, \f$ R^{T} \f$.
   *
   * `orbs` indexes positions within the \f$ (2l+1) \f$ shell, in the orbital ordering of the DFT code, and is
   * ordered: entry `i` names the shell orbital that is correlated orbital `i`. The selection need not be
   * contiguous -- for VASP's \f$ l=2 \f$ ordering (xy, yz, z², xz, x²-y²), t2g is `{0,1,3}` and \f$ e_g \f$ is
   * `{2,4}`.
   *
   * If the correlated space spans the whole shell (`dim == 2l+1`), `orbs` may be omitted. Otherwise it is
   * required: which orbitals the projectors span is not recorded by the DFT converters and cannot be inferred.
   *
   * \f$ R \f$ is read from `C_space`, i.e. it is the rotation the loader applied. `rotate_local_basis` does not
   * record itself there, so if you rotated the one-body elements yourself, rotate the tensor yourself too, with
   * @ref triqs::rotate_U_matrix_slater and the transpose of the rotation you passed.
   *
   * @param U_tensor Coulomb tensor spanning the full shell, e.g. from @ref triqs::slater_tensor.
   * @param C_space The local \f$ \mathcal{C} \f$ space.
   * @param atom Index of the atomic shell within `C_space`.
   * @param orbs Shell positions of the correlated orbitals; defaults to the whole shell.
   * @return Coulomb tensor of shape \f$ \mathrm{dim}^4 \f$ in the local basis.
   */
  nda::array<dcomplex, 4> to_local_basis(nda::array<dcomplex, 4> const &U_tensor, modest::local_space const &C_space, long atom,
                                         std::optional<std::vector<long>> const &orbs = {});

  ///@}

  /** @name Hamiltonians
   * Functions to create the impurity interaction using many-body operators.
   * @{
   */

  using op_name = std::pair<std::string, int>;

  operators::many_body_operator rename_op(operators::many_body_operator const &op, std::map<op_name, op_name> const &map);

  struct kanamori_params {
    bool spin_flip    = true;
    bool pair_hopping = true;
  };

  operators::many_body_operator h_int_kanamori(nda::matrix<double> const &U, nda::matrix<double> const &Uprime, double J_hund, long n_orb,
                                               std::vector<std::string> const &spin_names, kanamori_params const &params = {});

  operators::many_body_operator h_int_density(nda::matrix<double> const &U, nda::matrix<double> const &Uprime, double J_hund, long n_orb,
                                              std::vector<std::string> const &spin_names);

  operators::many_body_operator h_int_slater(nda::array<dcomplex, 4> const &U, long n_orb, std::vector<std::string> const &spin_names);

  /**
   * @ingroup hamiltonian
   * @brief Construct a density-density interaction Hamiltonian for one impurity.
   *
   * @details Create a density-density Hamiltonian
   * \f[
   *   H_{\mathrm{int}} = \frac{1}{2} \sum_{(i\sigma)\neq(j\sigma^{\prime})} U_{ij}^{\sigma\sigma^{\prime}}n_{i\sigma}
   *   n_{j\sigma^{\prime}}.
   * \f]
   *
   * Equivalent to @ref triqs::make_kanamori with the spin-flip and pair-hopping terms switched off.
   *
   * @param E The embedding.
   * @param imp_idx Index of the impurity.
   * @param U_int Hubbard \f$ U \f$.
   * @param U_prime \f$ U' \f$ (typically \f$ U' = U - 2J \f$).
   * @param J_hund Kanamori \f$ J \f$.
   * @return Many-body operator representing the Hamiltonian.
   */
  operators::many_body_operator make_density_density(modest::embedding const &E, long imp_idx, double U_int, double U_prime, double J_hund);

  /**
   * @ingroup hamiltonian
   * @brief Construct a Hubbard-Kanamori Hamiltonian for one impurity.
   *
   * @details Create a Hubbard-Kanamori Hamiltonian using the density-density, spin-flip, and pair-hopping interactions,
   * \f[
   *   H_{\mathrm{int}} = \frac{1}{2} \sum_{(i\sigma)\neq(j\sigma^{\prime})} U_{ij}^{\sigma\sigma^{\prime}}n_{i\sigma}
   *   n_{j\sigma^{\prime}} - \sum_{i\neq j}Jc_{i\uparrow}^{\dagger}c_{i\downarrow}c_{j\downarrow}^{\dagger}
   *   c_{j\uparrow} + \sum_{i\neq j} J c_{i\uparrow}^{\dagger}c_{i\downarrow}^{\dagger}c_{j\downarrow}c_{j\uparrow}.
   * \f]
   *
   * The operator names follow the block structure of impurity `imp_idx`, i.e. `(tau_gamma, orbital)`.
   *
   * @param E The embedding.
   * @param imp_idx Index of the impurity.
   * @param U_int  Hubbard \f$ U \f$.
   * @param U_prime \f$ U' \f$ (typically \f$ U' = U - 2J \f$).
   * @param J_hund Kanamori \f$ J \f$.
   * @param spin_flip Spin flip term.
   * @param pair_hopping Pair-hopping term.
   * @return Many-body operator representing the Hamiltonian.
   */
  operators::many_body_operator make_kanamori(modest::embedding const &E, long imp_idx, double U_int, double U_prime, double J_hund,
                                              bool spin_flip = true, bool pair_hopping = true);

  /**
   * @ingroup hamiltonian
   * @brief Construct a Slater Hamiltonian for one impurity from a Coulomb tensor.
   *
   * @details Create a Slater Hamiltonian using fully rotationally-invariant four-index interactions:
   * \f[
   *   H_{\mathrm{int}} = \frac{1}{2} \sum_{ijkl, \sigma\sigma^{\prime}} U_{ijkl}c^{\dagger}_{i\sigma}
   *   c^{\dagger}_{j\sigma^{\prime}}c_{l\sigma^{\prime}}c_{k\sigma}.
   * \f]
   *
   * `U_tensor` may be given at either of two sizes:
   *
   * * sized for impurity `imp_idx` itself, the usual case when the impurity is a whole atomic shell;
   * * sized for the whole correlated space \f$ \mathcal{C} \f$, in which case this impurity's block is cut out
   *   through the embedding. This is what @ref triqs::to_local_basis returns when the impurity is a part of
   *   \f$ \mathcal{C} \f$ split off with `split_imp`, e.g. a t2g impurity from a full d shell. It requires the
   *   impurity's orbitals to form one contiguous block of \f$ \mathcal{C} \f$; if they do not, restrict the
   *   tensor yourself and pass the impurity-sized one.
   *
   * @param E The embedding.
   * @param imp_idx Index of the impurity.
   * @param U_tensor Coulomb tensor in the local basis, of shape \f$ n_{\mathrm{orb}}^4 \f$ or
   * \f$ \mathrm{dim}(\mathcal{C})^4 \f$, e.g. built with @ref triqs::slater_tensor followed by
   * @ref triqs::to_local_basis.
   * @return Many-body operator representing the Hamiltonian.
   */
  operators::many_body_operator make_slater(modest::embedding const &E, long imp_idx, nda::array<dcomplex, 4> const &U_tensor);
  ///@}

} // namespace triqs
