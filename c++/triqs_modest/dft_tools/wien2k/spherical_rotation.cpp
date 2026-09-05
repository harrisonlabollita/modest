// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "../spherical_rotation.hpp"
#include <stdexcept>

namespace triqs::modest::dft_tools::wien2k {

  nda::matrix<dcomplex> get_spherical_to_dft_rotation(long l) {
    if (l == 0) throw std::runtime_error("No l=0 implementation");

    auto sqrt2 = std::numbers::sqrt2;
    auto I     = dcomplex(0, 1.0);
    auto size  = 2 * l + 1;
    auto Ylm   = nda::zeros<dcomplex>(size, size);
    // l = 1 : "x", "y", "z"
    if (l == 1) {
      Ylm(0, 0) = 1.0 / sqrt2;
      Ylm(1, 0) = I / sqrt2;
      Ylm(2, 1) = 1.0;
      Ylm(0, 2) = -1.0 / sqrt2;
      Ylm(1, 2) = I / sqrt2;
    }
    // l = 2: "z^2", "x^2-y^2", "xy", "yz", "xz"
    else if (l == 2) {
      Ylm(0, 2) = 1.0;
      Ylm(1, 0) = 1.0 / sqrt2;
      Ylm(2, 0) = -1.0 / sqrt2;
      Ylm(3, 1) = 1.0 / sqrt2;
      Ylm(4, 1) = 1.0 / sqrt2;
      Ylm(1, 4) = 1.0 / sqrt2;
      Ylm(2, 4) = 1.0 / sqrt2;
      Ylm(3, 3) = -1.0 / sqrt2;
      Ylm(4, 3) = 1.0 / sqrt2;
    }
    // There is no single Y_lm -> Wien2k rotation for an f shell: dmftproj builds a point-group-specific
    // transformation per site (see dmftproj/SRC_templates/case.cf_f_mm2). Returning a generic matrix here would
    // silently give the wrong Coulomb tensor, so refuse instead. TRIQS' spherical_to_cubic does the same.
    else if (l == 3) {
      throw std::runtime_error(
         "[wien2k::get_spherical_to_dft_rotation] No Wien2k f-shell (l=3) rotation: dmftproj uses point-group-specific "
         "transformations per site (see dmftproj/SRC_templates/case.cf_f_mm2). Pass the rotation explicitly.");
    } else
      throw std::runtime_error(fmt::format("[wien2k::get_spherical_to_dft_rotation] No implementation for l = {}.", l));
    return Ylm;
  }
} // namespace triqs::modest::dft_tools::wien2k
