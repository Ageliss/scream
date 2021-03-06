#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>


namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocAdvSgsTke {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;

    // Tests for the subroutine adv_sgs_tke
    //   in the SHOC TKE module.

    // For this routine all inputs are on midpoint grid and
    //  there are no vertical derivatives, therefore we will
    //  consider one vertical level per test. Each column will
    //  be loaded up with a different test (growth, decay)

    // FIRST TEST
    // TKE growth/decay test.  Given the correct conditions, verify
    //  that TKE grows/decays over a timestep.  For this we choose a
    //  large/small mixing length, positive buoyancy flux, positive
    //  shear term, and large eddy diffusivity.

    // Values in the FIRST column represents the growth test
    // Values in the SECOND column represents the decay test

    // Define timestep [s]
    static constexpr Real dtime = 20.0;
    // SHOC mixing length [m]
    static constexpr Real shoc_mix_gr[shcol] = {20000.0, 20.};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec_gr[shcol] = {0.5, -0.5};
    // Shear production term [s-2]
    static constexpr Real sterm_gr[shcol] = {0.5, 0.0};
    // TKE initial value
    Real tke_init_gr[shcol] = {0.004, 0.4};

    // Initialize data structure for bridgeing to F90
    SHOCAdvsgstkeData SDS(shcol, nlev, dtime);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.dtime == dtime) );
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.shoc_mix[offset] = shoc_mix_gr[s];
	SDS.wthv_sec[offset] = wthv_sec_gr[s];
	SDS.sterm_zt[offset] = sterm_gr[s];
	SDS.tke[offset] = tke_init_gr[s];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;
	// time step, mixing length, TKE values,
	// shear terms should all be greater than zero
	REQUIRE(SDS.dtime > 0.0);
	REQUIRE(SDS.shoc_mix[offset] > 0.0);
	REQUIRE(SDS.tke[offset] > 0.0);
	REQUIRE(SDS.sterm_zt[offset] >= 0.0);
      }
    }

    // Call the fortran implementation
    adv_sgs_tke(SDS);

    // Check to make sure that there has been
    //  TKE growth
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	if (s == 0){
          // Growth check
          REQUIRE(SDS.tke[offset] > tke_init_gr[s]);
	}
	else{
          // Decay check
          REQUIRE(SDS.tke[offset] < tke_init_gr[s]);
	}
      }
    }

    // SECOND TEST
    // TKE Dissipation test.  Given input values that are identical
    //  in two columns, verify that the dissipation rate is higher
    //  when the length scale is lower.

    // SHOC mixing length [m]
    static constexpr Real shoc_mix_diss[shcol] = {1000.0, 500.};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec_diss = 0.1;
    // Shear production term [s-2]
    static constexpr Real sterm_diss = 0.01;
    // TKE initial value
    Real tke_init_diss= 0.1;

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.shoc_mix[offset] = shoc_mix_diss[s];
	SDS.wthv_sec[offset] = wthv_sec_diss;
	SDS.sterm_zt[offset] = sterm_diss;
	SDS.tke[offset] = tke_init_diss;
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;
	// time step, mixing length, TKE values,
	// shear terms should all be greater than zero
	REQUIRE(SDS.dtime > 0.0);
	REQUIRE(SDS.shoc_mix[offset] > 0.0);
	REQUIRE(SDS.tke[offset] > 0.0);
	REQUIRE(SDS.sterm_zt[offset] >= 0.0);
      }
    }

    // Call the fortran implementation
    adv_sgs_tke(SDS);

    // Check to make sure that the column with
    //  the smallest length scale has larger
    //  dissipation rate
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;
	// Get value corresponding to next column
	const auto offsets = n + (s+1) * nlev;
	if(SDS.shoc_mix[offset] > SDS.shoc_mix[offsets]){
          REQUIRE(SDS.a_diss[offset] < SDS.a_diss[offsets]);
	}
	else {
          REQUIRE(SDS.a_diss[offset] > SDS.a_diss[offsets]);
	}
      }
    }
  }

  static void run_bfb()
  {
    // TODO
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_adv_sgs_tke_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAdvSgsTke;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_adv_sgs_tke_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAdvSgsTke;

  TestStruct::run_bfb();
}

} // namespace
