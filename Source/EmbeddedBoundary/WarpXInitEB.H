/* Copyright 2021 Lorenzo Giacomel, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H
#define WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H

#include "EmbeddedBoundary/WarpXFaceInfoBox.H"

#include <AMReX_EBFabFactory.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>


#include <array>
#include <memory>

namespace warpx::embedded_boundary
{

    /**
    * \brief Compute the length of the mesh edges. Here the length is a value in [0, 1].
    *        An edge of length 0 is fully covered.
    */
    void ComputeEdgeLengths (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                                    const amrex::EBFArrayBoxFactory& eb_fact);
    /**
    * \brief Compute the area of the mesh faces. Here the area is a value in [0, 1].
    *        An edge of area 0 is fully covered.
    */
    void ComputeFaceAreas (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                                  const amrex::EBFArrayBoxFactory& eb_fact);
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


    class EmbeddedBoundary
    {

    public:

        /**
        * \brief Initialize information for cell extensions.
        *        The flags convention for m_flag_info_face is as follows
        *          - 0 for unstable cells
        *          - 1 for stable cells which have not been intruded
        *          - 2 for stable cells which have been intruded
        *        Here we cannot know if a cell is intruded or not so we initialize all stable cells with 1
        */
        void MarkCells (const std::array<amrex::Real,3>& cell_size, const int max_level);

        amrex::EBFArrayBoxFactory const& fieldEBFactory (int lev) const noexcept
        {
            return static_cast<amrex::EBFArrayBoxFactory const&>(*m_field_factory[lev]);
        }

    private:

        // Factory for field data
        amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_field_factory;

        //! EB: Lengths of the mesh edges
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_edge_lengths;
        //! EB: Areas of the mesh faces
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_face_areas;

        /** EB: for every mesh face flag_info_face contains a:
         *          * 0 if the face needs to be extended
         *          * 1 if the face is large enough to lend area to other faces
         *          * 2 if the face is actually intruded by other face
         * It is initialized in WarpX::MarkCells
         * This is only used for the ECT solver.*/
        amrex::Vector<std::array< std::unique_ptr<amrex::iMultiFab>, 3 > > m_flag_info_face;
        /** EB: for every mesh face face flag_ext_face contains a:
        *          * 1 if the face needs to be extended
        *          * 0 otherwise
        * It is initialized in WarpX::MarkCells and then modified in WarpX::ComputeOneWayExtensions
        * and in WarpX::ComputeEightWaysExtensions
        * This is only used for the ECT solver.*/
        amrex::Vector<std::array< std::unique_ptr<amrex::iMultiFab>, 3 > > m_flag_ext_face;

        /** EB: m_area_mod contains the modified areas of the mesh faces, i.e. if a face is enlarged it
        * contains the area of the enlarged face
        * This is only used for the ECT solver.*/
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_area_mod;
        /** EB: m_borrowing contains the info about the enlarged cells, i.e. for every enlarged cell it
        * contains the info of which neighbors are being intruded (and the amount of borrowed area).
        * This is only used for the ECT solver.*/
        amrex::Vector<std::array< std::unique_ptr<amrex::LayoutData<FaceInfoBox> >, 3 > > m_borrowing;

    };



}

#endif //WARPX_EMBEDDEDBOUNDARY_WARPXINITEB_H
