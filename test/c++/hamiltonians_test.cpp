#include <nda/gtest_tools.hpp>
#include <string>
#include "./common.hpp"
#include "triqs_modest/hamiltonians.hpp"
#include "triqs_modest/utils/h5_proxy.hpp"
#include "triqs_modest/dft_tools/utils.hpp"
#include "triqs_modest/dft_tools/spherical_rotation.hpp"

using namespace triqs::operators;

TEST(hamiltonian_tests, kanamori) { // NOLINT
  auto spin_names    = std::vector<std::string>{"up", "down"};
  auto n_orb         = 2;
  auto U             = 3.0;
  auto J             = 1.0;
  auto Up            = U - 2 * J;
  auto [Umat, Upmat] = U_matrix_kanamori(n_orb, U, Up, J);
  auto h_int         = h_int_kanamori(Umat, Upmat, J, n_orb, spin_names);

  auto h_ref = 3 * n("up", 0) * n("down", 0) + 3 * n("up", 1) * n("down", 1) + 3 * n("up", 2) * n("down", 2);

  // TODO: refactor reference data to explicit construction
  EXPECT_PRINT(
     "1*c_dag('down',0)*c_dag('up',0)*c('up',1)*c('down',1) + 3*c_dag('down',0)*c_dag('up',0)*c('up',0)*c('down',0) + 1*c_dag('down',0)*c_dag('up',1)*c('up',1)*c('down',0) + 1*c_dag('down',0)*c_dag('up',1)*c('up',0)*c('down',1) + 1*c_dag('down',1)*c_dag('up',0)*c('up',1)*c('down',0) + 1*c_dag('down',1)*c_dag('up',0)*c('up',0)*c('down',1) + 3*c_dag('down',1)*c_dag('up',1)*c('up',1)*c('down',1) + 1*c_dag('down',1)*c_dag('up',1)*c('up',0)*c('down',0)",
     h_int);
}

TEST(hamiltonian_tests, kanamori_with_rename) { // NOLINT
  auto spin_names    = std::vector<std::string>{"up", "down"};
  auto n_orb         = 2;
  auto U             = 3.0;
  auto J             = 1.0;
  auto Up            = U - 2 * J;
  auto [Umat, Upmat] = U_matrix_kanamori(n_orb, U, Up, J);

  std::map<op_name, op_name> op_map;
  for (auto const &s : spin_names) {
    for (auto o : range(n_orb)) { op_map[{s, o}] = {fmt::format("{}_{}", s, o), 0}; }
  }

  // TODO: refactor reference data to explicit construction
  auto h_int = rename_op(h_int_kanamori(Umat, Upmat, J, n_orb, spin_names), op_map);
  EXPECT_PRINT(
     "1*c_dag('down_0',0)*c_dag('up_0',0)*c('up_1',0)*c('down_1',0) + 3*c_dag('down_0',0)*c_dag('up_0',0)*c('up_0',0)*c('down_0',0) + 1*c_dag('down_0',0)*c_dag('up_1',0)*c('up_1',0)*c('down_0',0) + 1*c_dag('down_0',0)*c_dag('up_1',0)*c('up_0',0)*c('down_1',0) + 1*c_dag('down_1',0)*c_dag('up_0',0)*c('up_1',0)*c('down_0',0) + 1*c_dag('down_1',0)*c_dag('up_0',0)*c('up_0',0)*c('down_1',0) + 3*c_dag('down_1',0)*c_dag('up_1',0)*c('up_1',0)*c('down_1',0) + 1*c_dag('down_1',0)*c_dag('up_1',0)*c('up_0',0)*c('down_0',0)",
     h_int);
}

TEST(hamiltonian_test, kanamori_from_embedding) {

  auto spin_names    = std::vector<std::string>{"up", "down"};
  auto n_orb         = 3;
  auto U             = 3.0;
  auto J             = 0.5;
  auto Up            = U - 2 * J;
  auto [Umat, Upmat] = U_matrix_kanamori(n_orb, U, Up, J);

  std::map<op_name, op_name> op_map;
  for (auto const &s : spin_names) {
    for (auto o : range(n_orb)) { op_map[{s, o}] = {fmt::format("{}_{}", s, o), 0}; }
  }
  auto h_int1 = rename_op(h_int_kanamori(Umat, Upmat, J, n_orb, spin_names), op_map);

  std::string filename = "ref_data/svo-wien2k.ref.h5";
  auto [_, obe]        = one_body_elements_from_dft_converter(filename, 1.e-3);
  auto E               = make_embedding(obe.C_space);
  auto h_int2          = make_kanamori(E.sigma_names(), E.imp_decomposition(0), 3.0, 2.0, 0.5);

  EXPECT_TRUE((h_int1 - h_int2).is_zero());
}

// Kanamori averages of a Coulomb tensor: (U_avg, J_avg). These are basis independent for a full shell, so they
// must return U_int and J_hund whatever orbital basis the tensor was rotated into.
std::pair<double, double> kanamori_averages(nda::array<dcomplex, 4> const &U) {
  auto N       = U.extent(0);
  auto U_avg   = 0.0;
  auto UmJ_avg = 0.0;
  for (auto [m, mp] : product(range(N), range(N))) {
    U_avg += real(U(m, mp, m, mp));
    if (m != mp) UmJ_avg += real(U(m, mp, m, mp)) - real(U(m, mp, mp, m));
  }
  U_avg /= double(N * N);
  UmJ_avg /= double(N * (N - 1));
  return {U_avg, U_avg - UmJ_avg};
}

TEST(hamiltonian_tests, slater_tensor_uses_shell_l_not_orbital_count) {
  // SrVO3 t2g: the correlated space has 3 orbitals but the shell is l = 2. The tensor must be built for the
  // full 5-orbital shell -- inferring l from the orbital count would give l = 1 and the wrong physics.
  auto [_, obe]  = one_body_elements_from_dft_converter("ref_data/svo-wien2k.ref.h5", 1.e-3);
  auto const &C  = obe.C_space;
  auto const &sh = C.atomic_shells()[0];
  EXPECT_EQ(sh.dim, 3);
  EXPECT_EQ(sh.l, 2);

  auto U_shell = slater_tensor(C, 0, 3.0, 0.5);
  for (int ax = 0; ax < 4; ++ax) EXPECT_EQ(U_shell.extent(ax), 5);

  auto [U_avg, J_avg] = kanamori_averages(U_shell);
  EXPECT_NEAR(U_avg, 3.0, 1e-12);
  EXPECT_NEAR(J_avg, 0.5, 1e-12);

  // identical to building it by hand at l = 2 with the stored rotation
  auto ref = U_matrix_slater_local(2, C.rotation_from_spherical_to_dft_basis()(0), 3.0, 0.5);
  EXPECT_ARRAY_NEAR(U_shell, ref, 1e-14);
}

TEST(hamiltonian_tests, to_local_basis_requires_orbitals_for_a_subshell) {
  auto [_, obe] = one_body_elements_from_dft_converter("ref_data/svo-wien2k.ref.h5", 1.e-3);
  auto const &C = obe.C_space;
  auto U_shell  = slater_tensor(C, 0, 3.0, 0.5);

  // dim (3) < 2l+1 (5): the shell orbitals spanned by the projectors cannot be inferred
  EXPECT_THROW(to_local_basis(U_shell, C, 0), std::runtime_error);
  EXPECT_THROW(to_local_basis(U_shell, C, 0, std::vector<long>{0, 1}), std::runtime_error);    // wrong count
  EXPECT_THROW(to_local_basis(U_shell, C, 0, std::vector<long>{0, 1, 5}), std::runtime_error); // out of shell

  auto orbs  = std::vector<long>{0, 1, 3};
  auto U_loc = to_local_basis(U_shell, C, 0, orbs);
  for (int ax = 0; ax < 4; ++ax) EXPECT_EQ(U_loc.extent(ax), 3);

  // equivalent to selecting the orbitals and then applying transpose(R) by hand
  auto U_sel = nda::zeros<dcomplex>(3, 3, 3, 3);
  for (auto [i, j, k, m] : product(range(3), range(3), range(3), range(3))) U_sel(i, j, k, m) = U_shell(orbs[i], orbs[j], orbs[k], orbs[m]);
  EXPECT_ARRAY_NEAR(U_loc, rotate_U_matrix_slater(U_sel, nda::matrix<dcomplex>{transpose(C.rotation_from_dft_to_local_basis()(0, 0))}), 1e-14);
}

TEST(hamiltonian_tests, slater_tensor_full_shell_for_several_l) {
  // NiO: a full d shell (l = 2) and a full p shell (l = 1). Both span their whole shell, so no orbital
  // selection is needed and both must return the input U and J.
  auto [_, obe] = one_body_elements_from_dft_converter("ref_data/nio.ref.h5", 1.e-3);
  auto const &C = obe.C_space;
  for (auto [atom, n_shell] : std::vector<std::pair<long, long>>{{0, 5}, {1, 3}}) {
    auto U_loc = to_local_basis(slater_tensor(C, atom, 3.0, 0.5), C, atom);
    EXPECT_EQ(U_loc.extent(0), n_shell);
    auto [U_avg, J_avg] = kanamori_averages(U_loc);
    EXPECT_NEAR(U_avg, 3.0, 1e-12);
    EXPECT_NEAR(J_avg, 0.5, 1e-12);
  }
}

TEST(hamiltonian_tests, radial_integrals_reject_unsupported_l) {
  // l = 0 and l > 3 used to return an all-zero tensor rather than complaining.
  EXPECT_THROW(U_matrix_slater_spherical(0, 3.0, 0.5), std::runtime_error);
  EXPECT_THROW(U_matrix_slater_spherical(4, 3.0, 0.5), std::runtime_error);
}

TEST(hamiltonian_tests, spherical_umatrix_slater_construction) {
  // auto spherical_to_cubic_conventions = {"wien2k", "qe", "vasp", "wannier90"};
  auto root   = h5::proxy{"ref_data/u_matrix_slater.ref.h5", 'r'};
  auto U_int  = as<double>(root["U_int"]);
  auto J_hund = as<double>(root["J_hund"]);
  auto ls     = std::vector<long>{1, 2, 3};
  for (auto l : ls) {
    auto Umat_ref = as<nda::array<double, 4>>(root[std::to_string(l)]["spherical"]);
    auto Umat     = U_matrix_slater_spherical(l, U_int, J_hund);
    EXPECT_ARRAY_NEAR(Umat, Umat_ref, 1e-14);
  }
}

TEST(hamiltonian_tests, umatrix_slater_construction_all_dft_conventions) {
  // Reference tensors come from TRIQS' spherical_to_cubic(l, convention) -- see ref_data/umatrix.py.
  auto root   = h5::proxy{"ref_data/u_matrix_slater.ref.h5", 'r'};
  auto U_int  = as<double>(root["U_int"]);
  auto J_hund = as<double>(root["J_hund"]);

  auto conventions = std::vector<std::pair<std::string, DFTCode>>{
     {"wien2k", DFTCode::Wien2k}, {"vasp", DFTCode::VASP}, {"qe", DFTCode::QuantumEspresso}, {"wannier90", DFTCode::W90}};

  for (auto l : std::vector<long>{1, 2, 3}) {
    for (auto const &[name, code] : conventions) {
      // Wien2k has no generic f-shell rotation; it must refuse rather than return a wrong one.
      if (l == 3 and code == DFTCode::Wien2k) {
        EXPECT_THROW(dft_tools::get_spherical_to_dft_rotation(code, l), std::runtime_error);
        continue;
      }
      auto Umat_ref = as<nda::array<dcomplex, 4>>(root[std::to_string(l)][name]);
      auto T        = dft_tools::get_spherical_to_dft_rotation(code, l);
      auto Umat     = U_matrix_slater_local(l, T, U_int, J_hund);
      EXPECT_ARRAY_NEAR(Umat, Umat_ref, 1e-14) << "l = " << l << ", convention = " << name;
    }
  }
}
