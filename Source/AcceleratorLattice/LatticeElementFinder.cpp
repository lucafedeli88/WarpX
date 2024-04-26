/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "LatticeElementFinder.H"

#include "AcceleratorLattice.H"
#include "LatticeElements/HardEdgedQuadrupole.H"
#include "LatticeElements/HardEdgedPlasmaLens.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

void
LatticeElementFinder::InitElementFinder (const int nz, const amrex::Real dz, const amrex::Real zmin,
                                         const amrex::Real time,
                                         const amrex::Real gamma_boost,
                                         AcceleratorLattice const& accelerator_lattice)
{
    m_nz = nz;
    m_dz = dz;
    m_gamma_boost = gamma_boost;
    m_uz_boost = std::sqrt(gamma_boost * gamma_boost - 1._prt)*PhysConst::c;

    AllocateIndices(accelerator_lattice);

    UpdateIndices(zmin, time, accelerator_lattice);
}

void
LatticeElementFinder::AllocateIndices (AcceleratorLattice const& accelerator_lattice)
{
    // Allocate the space for the indices for each element type.
    // Note that this uses m_nz since the information is saved per node.

    if (accelerator_lattice.h_quad.nelements > 0) {
        d_quad_indices.resize(m_nz);
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        d_plasmalens_indices.resize(m_nz);
    }
}

void
LatticeElementFinder::UpdateIndices (const amrex::Real zmin, const amrex::Real time,
                                     AcceleratorLattice const& accelerator_lattice)
{
    // Update the location of the index grid.
    m_zmin = zmin;
    m_time = time;

    if (accelerator_lattice.h_quad.nelements > 0) {
        setup_lattice_indices(accelerator_lattice.h_quad.d_zs,
                              accelerator_lattice.h_quad.d_ze,
                              d_quad_indices);
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        setup_lattice_indices(accelerator_lattice.h_plasmalens.d_zs,
                              accelerator_lattice.h_plasmalens.d_ze,
                              d_plasmalens_indices);
    }
}

LatticeElementFinderDevice
LatticeElementFinder::GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset,
                                              AcceleratorLattice const& accelerator_lattice) const
{
    LatticeElementFinderDevice result;
    result.InitLatticeElementFinderDevice(a_pti, a_offset, accelerator_lattice, *this);
    return result;
}

void
LatticeElementFinderDevice::InitLatticeElementFinderDevice (WarpXParIter const& a_pti, int const a_offset,
                                                            const amrex::Real gamma_boost, const amrex::Vector<amrex::Real>& dt_array,
                                                            AcceleratorLattice const& accelerator_lattice,
                                                            LatticeElementFinder const & h_finder)
{
    int const lev = a_pti.GetLevel();

    m_get_position = GetParticlePosition<PIdx>(a_pti, a_offset);
    const auto& attribs = a_pti.GetAttribs();
    m_ux = attribs[PIdx::ux].dataPtr() + a_offset;
    m_uy = attribs[PIdx::uy].dataPtr() + a_offset;
    m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
    m_dt = dt_array[lev];

    m_gamma_boost = gamma_boost;
    m_uz_boost = std::sqrt(gamma_boost*gamma_boost - 1._prt)*PhysConst::c;

    m_zmin = h_finder.m_zmin;
    m_dz = h_finder.m_dz;
    m_time = h_finder.m_time;

    if (accelerator_lattice.h_quad.nelements > 0) {
        d_quad = accelerator_lattice.h_quad.GetDeviceInstance();
        d_quad_indices_arr = h_finder.d_quad_indices.data();
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        d_plasmalens = accelerator_lattice.h_plasmalens.GetDeviceInstance();
        d_plasmalens_indices_arr = h_finder.d_plasmalens_indices.data();
    }

}

void
LatticeElementFinder::setup_lattice_indices (amrex::Gpu::DeviceVector<amrex::ParticleReal> const & zs,
                       amrex::Gpu::DeviceVector<amrex::ParticleReal> const & ze,
                       amrex::Gpu::DeviceVector<int> & indices) const
{

    using namespace amrex::literals;

    const auto nelements = static_cast<int>(zs.size());
    amrex::ParticleReal const * zs_arr = zs.data();
    amrex::ParticleReal const * ze_arr = ze.data();
    int * indices_arr = indices.data();

    amrex::Real const zmin = m_zmin;
    amrex::Real const dz = m_dz;

    amrex::ParticleReal const gamma_boost = m_gamma_boost;
    amrex::ParticleReal const uz_boost = m_uz_boost;
    amrex::Real const time = m_time;

    amrex::ParallelFor( m_nz,
        [=] AMREX_GPU_DEVICE (int iz) {

            // Get the location of the grid node
            amrex::Real z_node = zmin + iz*dz;

            if (gamma_boost > 1._prt) {
                // Transform to lab frame
                z_node = gamma_boost*z_node + uz_boost*time;
            }

            // Find the index to the element that is closest to the grid cell.
            // For now, this assumes that there is no overlap among elements of the same type.
            for (int ie = 0 ; ie < nelements ; ie++) {
                // Find the mid points between element ie and the ones before and after it.
                // The first and last element need special handling.
                const amrex::ParticleReal zcenter_left = (ie == 0)?
                    (std::numeric_limits<amrex::ParticleReal>::lowest()) : (0.5_prt*(ze_arr[ie-1] + zs_arr[ie]));
                const amrex::ParticleReal zcenter_right = (ie < nelements - 1)?
                    (0.5_prt*(ze_arr[ie] + zs_arr[ie+1])) : (std::numeric_limits<amrex::ParticleReal>::max());
                if (zcenter_left <= z_node && z_node < zcenter_right) {
                    indices_arr[iz] = ie;
                }
            }
        });
}
