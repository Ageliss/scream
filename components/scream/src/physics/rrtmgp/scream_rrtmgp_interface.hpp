#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "mo_gas_optics_rrtmgp.h"

namespace scream {
    namespace rrtmgp {
        /* 
         * Objects containing k-distribution information need to be initialized
         * once and then persist throughout the life of the program, so we
         * declare them here within the rrtmgp namespace.
         */
        extern GasOpticsRRTMGP k_dist_sw;
        extern GasOpticsRRTMGP k_dist_lw;
        extern void rrtmgp_initialize();
        extern void rrtmgp_main(
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, real2d &col_dry,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0);
        extern void rrtmgp_finalize();
        extern void rrtmgp_sw(
                GasOpticsRRTMGP &k_dist, 
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, real2d &col_dry, 
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0);
    } // namespace rrtmgp
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
