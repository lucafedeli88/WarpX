/* Copyright 2025      Andrew Myers, Ann Almgren, Aurore Blelly
 *                     Axel Huebl, Burlen Loring, David Grote
 *                     Glenn Richardson, Junmin Gu, Luca Fedeli
 *                     Mathieu Lobet, Maxence Thevenet, Michael Rowan
 *                     Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *                     Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_READ_PARAMETERS_H_
#define WARPX_READ_PARAMETERS_H_

#include <AMReX_REAL.H>

#include <limits>
#include <string>
#include <vector>

namespace warpx::initialization
{
    [[nodiscard]] std::string read_authors ();

    struct stop_condition {int max_step; amrex::Real stop_time;};
    [[nodiscard]] stop_condition read_stop_condition ();

    [[nodiscard]] std::string read_restart ();

    [[nodiscard]] bool has_checkpoint_diags ();

    namespace defaults
    {
        const auto default_max_step = std::numeric_limits<int>::max();
        const auto default_stop_time = std::numeric_limits<amrex::Real>::max();
    }
}

#endif //WARPX_READ_PARAMETERS_H_
