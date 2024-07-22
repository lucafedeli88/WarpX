/* Copyright 2021 Lorenzo Giacomel, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H
#define WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
#include <memory>

namespace warpx::embedded_boundary
{

    /**
    * \brief Scale the edges lengths by the mesh width to obtain the real lengths.
    */
    void ScaleEdges (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                            const std::array<amrex::Real,3>& cell_size);
    /**
    * \brief Scale the edges areas by the mesh width to obtain the real areas.
    */
    void ScaleAreas (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                            const std::array<amrex::Real,3>& cell_size);

}

#endif //WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H
