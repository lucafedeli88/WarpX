
/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PARTICLES_RADIATION_H
#define WARPX_PARTICLES_RADIATION_H

#include "RadiationHandler_fwd.H"

#include "Utils/Parser/IntervalsParser.H"

#include "Particles/WarpXParticleContainer_fwd.H"

#include <ablastr/math/Complex.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Scan.H>
#include <AMReX_Vector.H>

#include <memory>
#include <optional>
#include <string>

/* \brief CollisionHandler is a class that contains the
 * calculation of radiation for particles at each frequencies and angles.
 */
class RadiationHandler
{
public:
    RadiationHandler (const amrex::Array<amrex::Real,3>& center, const amrex::Geometry& geom, int shape_factor);

    /* Perform the calculation of the radiation */
    void add_radiation_contribution (
        const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc,
        amrex::Real current_time, int timestep);

    void dump_radiation (amrex::Real dt, int timestep, const std::string& filename);
    void Integral_overtime (const amrex::Real dt);

    enum struct Type
    {
        cartesian,
        spherical
    };

    void prepare_data_on_gpu ();

private:

    void gather_and_write_radiation (const std::string& filename, const int timestep);

    // Frequency range
    amrex::Array<amrex::Real,2> m_omega_range;
    int m_omega_points;
    amrex::Gpu::DeviceVector<amrex::Real> m_omegas;

    // Auxiliary coefficients for form factors
    amrex::Gpu::DeviceVector<amrex::Real> m_FF;

    // Dimensions of the detector
    amrex::Array<amrex::Real,2> m_ang_range;

    amrex::Real m_det_distance;
    amrex::Array<int,2> m_det_pts;
    amrex::Array<amrex::Real,3> m_center;
    amrex::Array<amrex::Real,3> m_det_direction;
    amrex::Array<amrex::Real,3> m_det_orientation;

    amrex::Array<amrex::Real,3> m_d;
    int m_shape_factor;

    amrex::Array<amrex::Vector<amrex::Real>,2> m_grid;

    amrex::Gpu::DeviceVector<amrex::Real> m_det_x;
    amrex::Gpu::DeviceVector<amrex::Real> m_det_y;
    amrex::Gpu::DeviceVector<amrex::Real> m_det_z;

    amrex::Gpu::DeviceVector<amrex::Real> m_det_n_x;
    amrex::Gpu::DeviceVector<amrex::Real> m_det_n_y;
    amrex::Gpu::DeviceVector<amrex::Real> m_det_n_z;

    amrex::Gpu::DeviceVector<ablastr::math::Complex> m_radiation_data;
    amrex::Gpu::DeviceVector<amrex::Real> m_radiation_calculation;

    Type m_radiation_type;

    utils::parser::IntervalsParser m_output_intervals_parser;

    bool m_use_logspace_for_omegas;

    std::optional<int> m_step_start = std::nullopt;
    std::optional<int> m_step_stop = std::nullopt;
    std::optional<int> m_step_skip = std::nullopt;

    std::optional<std::array<amrex::Real,2>> m_gamma_range = std::nullopt;
    amrex::Gpu::DeviceVector<int> m_mask;
    amrex::Gpu::DeviceVector<int> m_offset;
    amrex::Gpu::DeviceVector<int> m_idx;
};
#endif // WARPX_PARTICLES_RADIATION_H
