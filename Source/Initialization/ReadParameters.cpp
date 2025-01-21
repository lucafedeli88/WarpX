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

#include "ReadParameters.H"

#include "Utils/Parser/ParserUtils.H"

#include <AMReX_ParmParse.H>

#include <algorithm>

using namespace warpx::initialization;
using namespace amrex;

std::string read_authors ()
{
    const ParmParse pp;
    std::string authors;
    pp.query("authors", authors);
    return authors;
}

stop_condition read_stop_condition ()
{
    const ParmParse pp;
    int max_step = defaults::default_max_step;
    amrex::Real stop_time = defaults::default_stop_time;
    utils::parser::queryWithParser(pp, "max_step", max_step);
    utils::parser::queryWithParser(pp, "stop_time", stop_time);

    return {max_step, stop_time};
}

std::string read_restart ()
{
    const ParmParse pp_amr("amr");
    std::string restart_chkfile;
    pp_amr.query("restart", restart_chkfile);
    return restart_chkfile;
}

bool has_checkpoint_diags ()
{
    const ParmParse pp("diagnostics");
    std::vector<std::string> diags_names;
    pp.queryarr("diags_names", diags_names);

    const auto HasCheckPoint = [](const auto& name){
        std::string format;
        ParmParse{name}.query("format", format);
        return (format == "checkpoint");
      };

    return std::any_of(
        diags_names.cbegin(),
        diags_names.cend(),
        HasCheckPoint);
}
