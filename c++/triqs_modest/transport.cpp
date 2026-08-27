// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./transport.hpp"
#include "./lattice_gf_helpers.hpp"
#include "utils/defs.hpp"
#include <h5/h5.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

namespace triqs::modest {

  namespace detail {

    // Map a Cartesian direction character to an index.
    long dir_char_to_int(char c) {
      switch (c) {
        case 'x': return 0;
        case 'y': return 1;
        case 'z': return 2;
        default: throw std::runtime_error{fmt::format("transport: invalid direction character '{}'", c)};
      }
    }

    // Parse direction labels like "xy" into (alpha, beta) index pairs.
    std::vector<std::array<long, 2>> parse_directions(std::vector<std::string> const &directions) {
      std::vector<std::array<long, 2>> pairs;
      pairs.reserve(directions.size());
      for (auto const &d : directions) {
        if (d.size() != 2) throw std::runtime_error{fmt::format("transport: direction '{}' must be two characters (e.g. \"xy\")", d)};
        pairs.push_back({dir_char_to_int(d[0]), dir_char_to_int(d[1])});
      }
      return pairs;
    }

    // ω comes from mesh::refreq (uniform), so ω_n + Ω is a constant index shift. Ω is repinned to the
    // nearest multiple of dω and the repinned mesh is returned, so callers integrate against the Ω used.
    std::pair<nda::array<long, 1>, nda::array<double, 1>> build_omega_shifts(nda::array<double, 1> const &omega, nda::array<double, 1> const &Om) {
      long n_w = omega.extent(0);
      if (n_w < 2) throw std::runtime_error{"transport: the ω mesh needs at least two points to define a Ω shift."};
      double d_omega = (omega(n_w - 1) - omega(0)) / double(n_w - 1);

      long nq     = Om.extent(0);
      auto iOm    = nda::array<long, 1>(nq);
      auto Om_eff = nda::array<double, 1>(nq);
      for (long iq = 0; iq < nq; ++iq) {
        iOm(iq)    = std::llround(Om(iq) / d_omega);
        Om_eff(iq) = double(iOm(iq)) * d_omega;
      }
      return {iOm, Om_eff};
    }

    // The Cartesian symmetry operations to average over (fall back to the identity if none are provided).
    std::vector<nda::matrix<double>> symmetry_ops(band_velocities const &vel) {
      if (vel.rot_symmetries.empty()) return {nda::eye<double>(3)};
      return vel.rot_symmetries;
    }

    // Velocity block matrices v_α(k) for each Cartesian direction (already stored on the intersection window).
    std::array<nda::matrix<dcomplex>, 3> velocity_blocks(band_velocities const &vel, long sigma, long k_idx) {
      auto vfull = vel.v(sigma, k_idx); // (3, N_ν_v, N_ν_v)
      std::array<nda::matrix<dcomplex>, 3> v_block;
      for (long d = 0; d < 3; ++d) v_block[d] = nda::matrix<dcomplex>{vfull(d, r_all, r_all)};
      return v_block;
    }

    // Rotate the Cartesian velocity blocks by R:  vR[d] = Σ_c R(d,c) v_block[c].
    std::array<nda::matrix<dcomplex>, 3> rotate_velocity_blocks(std::array<nda::matrix<dcomplex>, 3> const &v_block, nda::matrix<double> const &R) {
      long n = v_block[0].extent(0);
      std::array<nda::matrix<dcomplex>, 3> vR;
      for (long d = 0; d < 3; ++d) {
        vR[d] = nda::zeros<dcomplex>(n, n);
        for (long c = 0; c < 3; ++c) vR[d] += R(d, c) * v_block[c];
      }
      return vR;
    }

  } // namespace detail

  transport_distribution_t transport_distribution(one_body_elements_on_grid const &obe, double mu,
                                                  block2_gf<mesh::refreq, matrix_valued> const &Sigma_w, nda::array<double, 1> const &Om_mesh,
                                                  std::vector<std::string> const &directions, double broadening) {
    if (!obe.velocities) throw std::runtime_error{"transport_distribution: obe.velocities is empty. Load with read_velocities=true."};
    if (obe.H.matrix_valued)
      throw std::runtime_error{"transport_distribution: matrix-valued H(k) is not supported (expected a diagonal band-basis dispersion)."};

    auto const &vel = *obe.velocities;
    if (vel.n_directions() != 3)
      throw std::runtime_error{fmt::format("transport_distribution: expected 3 Cartesian velocity directions, got {}", vel.n_directions())};
    double V       = *obe.cell_volume;
    auto spin_kind = obe.C_space.spin_kind();
    long n_spin    = (spin_kind == spin_kind_e::Polarized) ? 2 : 1;
    long spin_pol  = (spin_kind == spin_kind_e::NonPolarized) ? 0 : 1;

    auto const &mesh = Sigma_w(0, 0).mesh();
    long n_w         = mesh.size();
    auto im          = dcomplex(0, 1.0);
    auto delta       = im * broadening;

    // Real ω grid (from Σ's mesh) and the complex arguments ω + iδ used to build G.
    auto omega_arr = nda::array<double, 1>(n_w);
    for (auto &&[i, w] : enumerate(mesh)) omega_arr(i) = dcomplex(w).real();
    auto omegas = mesh | stdv::transform([&](auto w) { return dcomplex(w) + delta; }) | tl::to<std::vector>();

    auto dir_pairs     = detail::parse_directions(directions);
    long n_dir         = long(dir_pairs.size());
    long n_Om          = Om_mesh.extent(0);
    auto [iOm, Om_eff] = detail::build_omega_shifts(omega_arr, Om_mesh);

    // Cartesian directions actually used by the requested pairs (only these get a batched v·A GEMM below).
    std::array<bool, 3> need_dir{false, false, false};
    for (auto const &p : dir_pairs) {
      need_dir[p[0]] = true;
      need_dir[p[1]] = true;
    }

    // Active subspace and Σ_active per spin (independent of k).
    auto decomp       = get_struct(Sigma_w).dims(r_all, 0) | tl::to<std::vector>();
    auto active       = detail::detect_active_subspace(Sigma_w, decomp);
    auto Sa_per_sigma = range(n_spin) | stdv::transform([&](long s) { return detail::compute_sigma_active(Sigma_w, active, s); }) //
       | tl::to<std::vector>();

    auto syms       = detail::symmetry_ops(vel);
    double sym_norm = 1.0 / double(syms.size());

    auto Gamma       = nda::zeros<double>(n_dir, n_Om, n_w);
    auto Gamma_intra = nda::zeros<double>(n_dir, n_Om, n_w);
    auto Gamma_inter = nda::zeros<double>(n_dir, n_Om, n_w);

    mpi::communicator comm = {};
    long n_k               = obe.H.n_k();

    for (auto k_idx : mpi::chunk(range(n_k), comm)) {
      double w_k = obe.H.k_weights(k_idx);
      for (auto sigma : range(n_spin)) {
        // The velocities are stored on the intersection of the dispersion/A window and the velocity window, so
        // only its offset into the A array is needed here (see band_velocities).
        long a_off = vel.A_offset(sigma, k_idx);
        long n_ov  = vel.N_nu_v(sigma, k_idx);
        if (n_ov <= 0) continue;
        auto A_slice = nda::range(a_off, a_off + n_ov);

        // Full band-basis G, then the block spectral function A = i(G − G†)/(2π) on the intersection.
        auto G_band = detail::build_G_band(obe, active, Sa_per_sigma[sigma], omegas, mu, k_idx, sigma);
        // Spectral function A = i(G − G†)/(2π) on the intersection, stacked over ω as block-columns:
        //   A2(j, n·n_ov + l) = A_k(ω_n)_{jl}   (n_ov, n_ω·n_ov).
        // Stacking turns each v·A (a matrix product for every ω) into a single GEMM over all ω at once.
        auto A2 = nda::matrix<dcomplex>(n_ov, n_w * n_ov);
        for (auto n : range(n_w)) {
          auto Gsub                                       = nda::matrix<dcomplex>{G_band(n, A_slice, A_slice)};
          A2(r_all, nda::range(n * n_ov, (n + 1) * n_ov)) = nda::matrix<dcomplex>{im * (Gsub - dagger(Gsub)) / (2.0 * M_PI)};
        }

        auto v_block = detail::velocity_blocks(vel, sigma, k_idx); // per-direction (n_ov × n_ov)

        for (auto const &R : syms) {
          auto vR = detail::rotate_velocity_blocks(v_block, R); // vR[d] = Σ_c R(d,c) v_block[c]

          // Compute the velocity-weighted spectral functions for each direction.
          std::array<nda::array<dcomplex, 3>, 3> vA, vAd;
          for (auto d : range(3)) {
            if (!need_dir[d]) continue;
            auto vA2 = nda::matrix<dcomplex>{vR[d] * A2}; // (n_ov, n_ω·n_ov) — single GEMM over all ω
            vA[d]    = nda::array<dcomplex, 3>(n_w, n_ov, n_ov);
            vAd[d]   = nda::array<dcomplex, 3>(n_w, n_ov, n_ov);
            for (auto n : range(n_w)) {
              vA[d](n, r_all, r_all) = vA2(r_all, nda::range(n * n_ov, (n + 1) * n_ov));
              for (auto i : range(n_ov))
                for (auto jj : range(n_ov)) vAd[d](n, i, jj) = vR[d](i, i) * A2(i, n * n_ov + jj);
            }
          }

          for (auto ip : range(n_dir)) {
            long a = dir_pairs[ip][0], b = dir_pairs[ip][1];
            auto const &vA_a = vA[a];  // v_α   · A, direction α = a
            auto const &vA_b = vA[b];  // v_β   · A, direction β = b
            auto const &vd_a = vAd[a]; // v_α^d · A
            auto const &vd_b = vAd[b]; // v_β^d · A
            for (auto n : range(n_w)) {
              for (auto iq : range(n_Om)) {
                long j = n + iOm(iq);
                if (j < 0 || j >= n_w) continue;
                // total = Tr[v_α A(ω+Ω) v_β A(ω)]
                // intra = Tr[v_α^d A(ω+Ω) v_β^d A(ω)]   (diagonal velocity matrix elements)
                // inter = Tr[v_α^o A(ω+Ω) v_β^o A(ω)]   (off-diagonal, v^o = v − v^d)
                dcomplex tot = 0, dd = 0, oo = 0;
                for (auto i : range(n_ov)) {
                  for (auto jj : range(n_ov)) {
                    dcomplex a_f = vA_a(j, i, jj), a_d = vd_a(j, i, jj);
                    dcomplex b_f = vA_b(n, jj, i), b_d = vd_b(n, jj, i);
                    tot += a_f * b_f;
                    dd += a_d * b_d;
                    oo += (a_f - a_d) * (b_f - b_d);
                  }
                }
                double f = w_k * sym_norm;
                Gamma(ip, iq, n) += f * tot.real();
                Gamma_intra(ip, iq, n) += f * dd.real();
                Gamma_inter(ip, iq, n) += f * oo.real();
              }
            }
          }
        }
      }
    }

    Gamma       = mpi::all_reduce(Gamma);
    Gamma_intra = mpi::all_reduce(Gamma_intra);
    Gamma_inter = mpi::all_reduce(Gamma_inter);

    Gamma /= V;
    Gamma_intra /= V;
    Gamma_inter /= V;

    return {.Gamma             = std::move(Gamma),
            .Gamma_intra       = std::move(Gamma_intra),
            .Gamma_inter       = std::move(Gamma_inter),
            .omega_mesh        = std::move(omega_arr),
            .Om_mesh           = std::move(Om_eff),
            .directions        = directions,
            .spin_polarization = spin_pol};
  }

  transport_function_t transport_function(one_body_elements_on_grid const &obe, nda::array<double, 1> const &omega,
                                          std::vector<std::string> const &directions, double broadening) {
    if (!obe.velocities) throw std::runtime_error{"transport_function: obe.velocities is empty. Load with read_velocities=true."};
    if (!obe.cell_volume) throw std::runtime_error{"transport_function: obe.cell_volume is empty. Load with read_velocities=true."};

    auto const &vel = *obe.velocities;
    if (vel.n_directions() != 3)
      throw std::runtime_error{fmt::format("transport_function: expected 3 Cartesian velocity directions, got {}", vel.n_directions())};
    double V       = *obe.cell_volume;
    auto spin_kind = obe.C_space.spin_kind();
    long n_spin    = (spin_kind == spin_kind_e::Polarized) ? 2 : 1;

    auto dir_pairs = detail::parse_directions(directions);
    long n_dir     = long(dir_pairs.size());
    long n_w       = omega.extent(0);
    double eta2    = broadening * broadening;
    double norm    = broadening / M_PI;

    auto syms       = detail::symmetry_ops(vel);
    double sym_norm = 1.0 / double(syms.size());

    auto Phi = nda::zeros<double>(n_dir, n_w);

    mpi::communicator comm = {};
    long n_k               = obe.H.n_k();

    for (auto k_idx : mpi::chunk(range(n_k), comm)) {
      double w_k = obe.H.k_weights(k_idx);
      for (auto sigma : range(n_spin)) {
        // The velocities are stored on the intersection of the dispersion/A window and the velocity window, so
        // only its offset into the H/band array is needed here (see band_velocities).
        long a_off = vel.A_offset(sigma, k_idx);
        long n_ov  = vel.N_nu_v(sigma, k_idx);
        if (n_ov <= 0) continue;

        auto Hk      = obe.H.H(sigma, k_idx);                      // diagonal band energies (assumes band basis)
        auto v_block = detail::velocity_blocks(vel, sigma, k_idx); // per-direction (n_ov × n_ov)

        for (auto const &R : syms) {
          auto vR = detail::rotate_velocity_blocks(v_block, R); // vR[d] = Σ_c R(d,c) v_block[c]
          for (auto ip : range(n_dir)) {
            long a = dir_pairs[ip][0], b = dir_pairs[ip][1];
            for (auto ib : range(n_ov)) {
              double eps  = (Hk(a_off + ib, a_off + ib)).real(); // ε_n(k), measured from the same reference as omega
              double vv   = (vR[a](ib, ib) * vR[b](ib, ib)).real();
              double pref = w_k * sym_norm * vv;
              for (auto iw : range(n_w)) {
                double d2 = omega(iw) - eps;
                Phi(ip, iw) += pref * norm / (d2 * d2 + eta2);
              }
            }
          }
        }
      }
    }

    Phi = mpi::all_reduce(Phi);
    Phi /= V;

    return {.Phi = std::move(Phi), .omega_mesh = omega, .directions = directions};
  }

  // ------------------------------------------------------------------
  // h5 + printing
  // ------------------------------------------------------------------

  void h5_read(h5::group g, std::string const &name, transport_distribution_t &x) {
    auto sg = g.open_group(name);
    assert_hdf5_format(sg, x);
    h5_read(sg, "Gamma", x.Gamma);
    h5_read(sg, "Gamma_intra", x.Gamma_intra);
    h5_read(sg, "Gamma_inter", x.Gamma_inter);
    h5_read(sg, "omega_mesh", x.omega_mesh);
    h5_read(sg, "Om_mesh", x.Om_mesh);
    h5_read(sg, "directions", x.directions);
    h5_read(sg, "spin_polarization", x.spin_polarization);
  }
  void h5_write(h5::group g, std::string const &name, transport_distribution_t const &x) {
    auto sg = g.create_group(name);
    write_hdf5_format(sg, x);
    h5_write(sg, "Gamma", x.Gamma);
    h5_write(sg, "Gamma_intra", x.Gamma_intra);
    h5_write(sg, "Gamma_inter", x.Gamma_inter);
    h5_write(sg, "omega_mesh", x.omega_mesh);
    h5_write(sg, "Om_mesh", x.Om_mesh);
    h5_write(sg, "directions", x.directions);
    h5_write(sg, "spin_polarization", x.spin_polarization);
  }

  void h5_read(h5::group g, std::string const &name, transport_function_t &x) {
    auto sg = g.open_group(name);
    assert_hdf5_format(sg, x);
    h5_read(sg, "Phi", x.Phi);
    h5_read(sg, "omega_mesh", x.omega_mesh);
    h5_read(sg, "directions", x.directions);
  }
  void h5_write(h5::group g, std::string const &name, transport_function_t const &x) {
    auto sg = g.create_group(name);
    write_hdf5_format(sg, x);
    h5_write(sg, "Phi", x.Phi);
    h5_write(sg, "omega_mesh", x.omega_mesh);
    h5_write(sg, "directions", x.directions);
  }

  std::ostream &operator<<(std::ostream &out, transport_distribution_t const &x) {
    out << "Transport distribution Γ_αβ(ω,Ω):\n";
    out << fmt::format("  directions        = {}\n", fmt::join(x.directions, ", "));
    out << fmt::format("  Gamma shape       = {}\n", x.Gamma.shape());
    out << fmt::format("  n_omega, n_Om     = {}, {}\n", x.omega_mesh.extent(0), x.Om_mesh.extent(0));
    out << fmt::format("  spin_polarization = {}\n", x.spin_polarization);
    return out;
  }

  std::ostream &operator<<(std::ostream &out, transport_function_t const &x) {
    out << "Transport function Φ_αβ(ω):\n";
    out << fmt::format("  directions = {}\n", fmt::join(x.directions, ", "));
    out << fmt::format("  Phi shape  = {}\n", x.Phi.shape());
    return out;
  }

} // namespace triqs::modest
