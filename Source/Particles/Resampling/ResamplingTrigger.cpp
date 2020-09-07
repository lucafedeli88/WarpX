/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ResamplingTrigger.H"
#include "WarpX.H"

ResamplingTrigger::ResamplingTrigger ()
{
    amrex::ParmParse pprt("resampling_trigger");

    std::string resampling_trigger_int_string = "0";
    pprt.query("intervals", resampling_trigger_int_string);
    m_resampling_intervals = IntervalsParser(resampling_trigger_int_string);

    pprt.query("max_avg_ppc", m_max_avg_ppc);
}

bool ResamplingTrigger::triggered (const int timestep, const amrex::Real global_numparts) const
{
    if (!m_initialized) {initialize_global_numcells();};

    const amrex::Real avg_ppc = global_numparts/m_global_numcells;
    return (m_resampling_intervals.contains(timestep) ||
            avg_ppc > m_max_avg_ppc);
}

void ResamplingTrigger::initialize_global_numcells () const
{
    auto & warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.maxLevel(); lev++)
    {
        m_global_numcells +=  warpx.boxArray(lev).numPts();
    }
    m_initialized = true;
}
