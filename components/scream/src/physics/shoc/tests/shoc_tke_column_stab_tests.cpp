#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

TEST_CASE("shoc_tke_column_stab", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 5;

  // Tests for the subroutine integ_column_stability
  //   in the SHOC TKE module.

  // FIRST TEST
  // Symmetric test.  Verify that given a profile
  //  of symetric inputs for brunt vaisalla and dz profile
  //  that is uniform with height that the column integrated
  //  value is 0.0.

  // Define height thickness on nlev grid [m]
  //   Do a uniform grid for symetric test
  Real dz_zt[nlev] = {30., 30., 30., 30., 30.};
  // Define Pressure [hPa] (later convereted to Pa)
  Real pres[nlev] = {850., 900., 925.0, 950.0, 1000.0};
  // Brunt Vaisalla frequency [/s]
  Real brunt_sym[nlev] = {-0.5, -0.25, 0.0, 0.25, 0.5};

  // Convert pres to Pa
  for (Int n = 0; n < nlev; ++n){
    pres[n] = pres[n]*100.0;
  }

  // Initialzie data structure for bridgeing to F90
  SHOCColstabData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.shcol > 0);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.dz_zt[offset] = dz_zt[n];
      SDS.pres[offset] = pres[n];
      SDS.brunt[offset] = brunt_sym[n];
    }
  }

  // Check that the inputs make sense
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlev;
      // Should be greater than zero
      REQUIRE(SDS.dz_zt[offset] > 0.0);
      // Make sure all pressure levels are in the
      //  lower troposphere for this test
      REQUIRE(SDS.pres[offset] > 80000.0);
    }
  }

  // Call the fortran implementation
  integ_column_stability(SDS);

  // Check test
  //  Verify that output is zero
  for(Int s = 0; s < shcol; ++s) {
    REQUIRE(SDS.brunt_int[s] == 0.0);
  }

  // SECOND TEST
  // For set of inputs where brunt is negative at all
  //   points, then verify that output is negative

  // Brunt Vaisalla frequency [/s]
  Real brunt_neg[nlev] = {-0.3, -0.4, -0.1, -10.0, -0.5};

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      SDS.brunt[offset] = brunt_neg[n];
    }
  }

  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlev;
      // All points should be less than zero
      REQUIRE(SDS.brunt[offset] < 0.0);
    }
  }

  // Call the fortran implementation
  integ_column_stability(SDS);

  // Check test
  //  Verify that output is negative
  for(Int s = 0; s < shcol; ++s) {
    REQUIRE(SDS.brunt_int[s] < 0.0);
  }

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream