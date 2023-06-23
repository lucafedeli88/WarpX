/* Copyright 2019-2020 Andrew Myers, Burlen Loring, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Utils/Parser/ParserUtils.H"
#include "TextMsg.H"
#include "WarpXAlgorithmSelection.H"
#include "WarpXConst.H"
#include "WarpXProfilerWrapper.H"
#include "WarpXUtil.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <limits>

using namespace amrex;

/* \brief Function that sets the value of MultiFab MF to zero for z between
 * zmin and zmax.
 */
void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin, amrex::Real zmax){
    WARPX_PROFILE("WarpXUtil::NullifyMF()");
    int const ncomp = mf.nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();
        // Get box lower and upper physical z bound, and dz
        const amrex::Real zmin_box = WarpX::LowerCorner(bx, lev, 0._rt)[2];
        const amrex::Real zmax_box = WarpX::UpperCorner(bx, lev, 0._rt)[2];
        const amrex::Real dz  = WarpX::CellSize(lev)[2];
        // Get box lower index in the z direction
#if defined(WARPX_DIM_3D)
        const int lo_ind = bx.loVect()[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        const int lo_ind = bx.loVect()[1];
#else
        const int lo_ind = bx.loVect()[0];
#endif
        // Check if box intersect with [zmin, zmax]
        if ( (zmax>zmin_box && zmin<=zmax_box) ){
            const Array4<Real> arr = mf[mfi].array();
            // Set field to 0 between zmin and zmax
            ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept{
#if defined(WARPX_DIM_3D)
                    const Real z_gridpoint = zmin_box+(k-lo_ind)*dz;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    const Real z_gridpoint = zmin_box+(j-lo_ind)*dz;
#else
                    const Real z_gridpoint = zmin_box+(i-lo_ind)*dz;
#endif
                    if ( (z_gridpoint >= zmin) && (z_gridpoint < zmax) ) {
                        arr(i,j,k,n) = 0.;
                    }
                }
            );
        }
    }
}

namespace WarpXUtilIO{
    bool WriteBinaryDataOnFile(std::string filename, const amrex::Vector<char>& data)
    {
        std::ofstream of{filename, std::ios::binary};
        of.write(data.data(), data.size());
        of.close();
        return  of.good();
    }
}


namespace WarpXUtilLoadBalance
{
    bool doCosts (const amrex::LayoutData<amrex::Real>* costs, const amrex::BoxArray ba,
                  const amrex::DistributionMapping& dm)
    {
        const bool consistent = costs && (dm == costs->DistributionMap()) &&
            (ba.CellEqual(costs->boxArray())) &&
            (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers);
        return consistent;
    }
}
