/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "RadiationHandler.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"

#include "ablastr/warn_manager/WarnManager.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#ifdef WARPX_USE_OPENPMD
#   include "Diagnostics/OpenPMDHelpFunction.H"
#   include <openPMD/openPMD.hpp>
#endif

#include <cmath>
#include <tuple>
#include <vector>

using namespace amrex;
using namespace ablastr::math;
using namespace utils::parser;

#ifdef WARPX_USE_OPENPMD
    namespace io = openPMD;
#endif

auto const radiation_type_map  = std::map<std::string, RadiationHandler::Type>{
    {"cartesian", RadiationHandler::Type::cartesian},
    {"spherical", RadiationHandler::Type::spherical}
};

namespace
{

    auto compute_detector_positions(
        const amrex::Array<amrex::Real,3>& center,
        const amrex::Array<amrex::Real,3>& direction,
        const amrex::Real distance,
        const amrex::Array<amrex::Real,3>& orientation,
        const amrex::Array<int,2>& det_points,
        const amrex::Array<amrex::Real,2>& ang_range,
        const RadiationHandler::Type radiation_type)
    {

        WARPX_PROFILE("compute_detector_positions");

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            direction[0]*orientation[0] +
            direction[1]*orientation[1] +
            direction[2]*orientation[2] == 0,
            "Radiation detector orientation cannot be aligned with detector direction");

        const auto one_over_direction = 1.0_rt/std::sqrt(
            direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        const auto norm_direction = amrex::Array<amrex::Real,3>{
            direction[0]*one_over_direction,
            direction[1]*one_over_direction,
            direction[2]*one_over_direction};

        auto u = amrex::Array<amrex::Real,3>{
            orientation[0] - direction[0]*orientation[0],
            orientation[1] - direction[1]*orientation[1],
            orientation[2] - direction[2]*orientation[2]};
        const auto one_over_u = 1.0_rt/std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        u[0] *= one_over_u;
        u[1] *= one_over_u;
        u[2] *= one_over_u;

        auto v = amrex::Array<amrex::Real,3>{
            direction[1]*u[2]-direction[2]*u[1],
            direction[2]*u[0]-direction[0]*u[2],
            direction[0]*u[1]-direction[1]*u[0]};
        const auto one_over_v = 1.0_rt/std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        v[0] *= one_over_v;
        v[1] *= one_over_v;
        v[2] *= one_over_v;

        const auto how_many = det_points[0]*det_points[1];

        auto host_det_x = amrex::Vector<amrex::Real>(how_many);
        auto host_det_y = amrex::Vector<amrex::Real>(how_many);
        auto host_det_z = amrex::Vector<amrex::Real>(how_many);

        amrex::Array<amrex::Vector<amrex::Real>,2> grid;

        if(radiation_type == RadiationHandler::Type::cartesian){
#if defined(WARPX_DIM_3D)
            auto us = amrex::Vector<amrex::Real>(det_points[0]);
            const auto ulim = distance*std::tan(ang_range[0]*0.5_rt);
            amrex::linspace(us.begin(), us.end(), -ulim, ulim);
            auto vs = amrex::Vector<amrex::Real>(det_points[1]);
            const auto vlim = distance*std::tan(ang_range[1]*0.5_rt);
            amrex::linspace(vs.begin(), vs.end(), -vlim, vlim);
            for (int i = 0; i < det_points[0]; ++i)
            {
                for (int j = 0; j < det_points[1]; ++j)
                {
                    auto x = distance * norm_direction[0] + us[i]*u[0] + vs[j]*v[0];
                    auto y = distance * norm_direction[1] + us[i]*u[1] + vs[j]*v[1];
                    auto z = distance * norm_direction[2] + us[i]*u[2] + vs[j]*v[2];
                    host_det_x[i*det_points[1] + j] = center[0] + x;
                    host_det_y[i*det_points[1] + j] = center[1] + y;
                    host_det_z[i*det_points[1] + j] = center[2] + z;
                }
            }
            grid[0] = us;
            grid[1] = vs;
#else
    WARPX_ABORT_WITH_MESSAGE("Using cartesian radiation detector is only possible in 3D.");
#endif
        }
        if(radiation_type == RadiationHandler::Type::spherical){
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ang_range[0] >= 0 && ang_range[0] <= 2*ablastr::constant::math::pi,
                "angle_aperture[0] (phi) must be between 0 and 2*pi" );
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ang_range[1] >= 0 && ang_range[1] <= ablastr::constant::math::pi,
                "angle_aperture[1] (theta) must be between 0 and pi" );

            auto phis = amrex::Vector<amrex::Real>(det_points[0]);
            amrex::linspace(phis.begin(), phis.end(),-ang_range[0]*0.5_rt, ang_range[0]*0.5_rt);

#if defined(WARPX_DIM_3D)
            auto thetas = amrex::Vector<amrex::Real>(det_points[1]);
            amrex::linspace(thetas.begin(), thetas.end(),-ang_range[1]*0.5_rt, ang_range[1]*0.5_rt);
#else
            const auto thetas = amrex::Vector<amrex::Real>{0.0_rt};
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                det_points[1] == 1,
                "det_points[1].size() must be between 1 in 2D geometry" );
#endif

        for (int i = 0; i < det_points[0]; ++i)
        {
            for (int j = 0; j < det_points[1]; ++j)
            {
                const auto [sin_phi, cos_phi] = amrex::Math::sincos(phis[i]);
                const auto [sin_theta, cos_theta] = amrex::Math::sincos(thetas[j]);

                const auto wx = norm_direction[0]*cos_theta*cos_phi + u[0]*sin_theta + v[0]*cos_theta*sin_phi;
                const auto wy = norm_direction[1]*cos_theta*cos_phi + u[1]*sin_theta + v[1]*cos_theta*sin_phi;
                const auto wz = norm_direction[2]*cos_theta*cos_phi + u[2]*sin_theta + v[2]*cos_theta*sin_phi;

                const int idx = i*det_points[1] + j;
                host_det_x[idx] = center[0] + distance*wx;
                host_det_y[idx] = center[1] + distance*wy;
                host_det_z[idx] = center[2] + distance*wz;
            }
        }

        grid[0] = phis;
        grid[1] = thetas;
    }

        return std::make_tuple(host_det_x, host_det_y, host_det_z, grid);
    }
}




RadiationHandler::RadiationHandler(const amrex::Array<amrex::Real,3>& center, const amrex::Geometry& geom, const int shape_factor)
{
    WARPX_PROFILE("RadiationHandler::RadiationHandler");

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D)
    WARPX_ABORT_WITH_MESSAGE("Radiation is not supported yet in RZ and 1D.");
#endif

    // Read in radiation input
    const amrex::ParmParse pp_radiation("radiation");

    //type of detector
    std::string radiation_type = "spherical";
    pp_radiation.query("detector_type", radiation_type);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        radiation_type_map.find(radiation_type) != radiation_type_map.end(),
        radiation_type + " detector type is not supported yet");
    m_radiation_type = radiation_type_map.at(radiation_type);

    //Resolution in frequency of the detector
    auto omega_range = std::vector<amrex::Real>(2);
    getArrWithParser(pp_radiation, "omega_range", omega_range);
    std::copy(omega_range.begin(), omega_range.end(), m_omega_range.begin());
    getWithParser(pp_radiation, "omega_points", m_omega_points);
    //Angle theta AND phi
    auto ang_range = std::vector<amrex::Real>(2);
    getArrWithParser(pp_radiation, "angle_aperture", ang_range);
    std::copy(ang_range.begin(), ang_range.end(), m_ang_range.begin());

    //Detector parameters
    auto det_pts = std::vector<int>(2);
    getArrWithParser(pp_radiation, "detector_number_points", det_pts);
    std::copy(det_pts.begin(), det_pts.end(), m_det_pts.begin());
#if defined(WARPX_DIM_XZ)
    m_det_pts[1] = 1;
#endif

    auto det_direction = std::vector<amrex::Real>(3);
    getArrWithParser(pp_radiation, "detector_direction", det_direction);
#if defined(WARPX_DIM_XZ)
    det_direction[1] = 0.0;
#endif
    std::copy(det_direction.begin(), det_direction.end(), m_det_direction.begin());

#if defined(WARPX_DIM_XZ)
    m_det_orientation = amrex::Array<amrex::Real,3>{0.0,1.0,0.0};
#else
    auto det_orientation = std::vector<amrex::Real>(3);
    getArrWithParser(pp_radiation, "detector_orientation", det_orientation);
    std::copy(det_orientation.begin(), det_orientation.end(), m_det_orientation.begin());
#endif

    getWithParser(pp_radiation, "detector_distance", m_det_distance);

    const auto [det_x, det_y, det_z, grid] = compute_detector_positions(
        center, m_det_direction, m_det_distance,
        m_det_orientation, m_det_pts, m_ang_range, m_radiation_type);
    m_grid = grid;

    m_det_x = amrex::Gpu::DeviceVector<amrex::Real>(det_x.size());
    m_det_y = amrex::Gpu::DeviceVector<amrex::Real>(det_y.size());
    m_det_z = amrex::Gpu::DeviceVector<amrex::Real>(det_z.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
        det_x.begin(), det_x.end(), m_det_x.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
        det_y.begin(), det_y.end(), m_det_y.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
        det_z.begin(), det_z.end(), m_det_z.begin());

    // Pre-compute detector directions
    m_det_n_x = amrex::Gpu::DeviceVector<amrex::Real>(det_x.size());
    m_det_n_y = amrex::Gpu::DeviceVector<amrex::Real>(det_y.size());
    m_det_n_z = amrex::Gpu::DeviceVector<amrex::Real>(det_z.size());

    constexpr auto ncomp = 3;
    m_radiation_data = amrex::Gpu::DeviceVector<ablastr::math::Complex>(m_det_pts[0]*m_det_pts[1]*m_omega_points*ncomp);

    int t_use_logspace_for_omegas = 0;
    pp_radiation.query("use_logspace_for_omegas", t_use_logspace_for_omegas);
    m_use_logspace_for_omegas = static_cast<bool>(t_use_logspace_for_omegas);

    auto t_omegas = amrex::Vector<amrex::Real>(m_omega_points);
    if (m_use_logspace_for_omegas){
        amrex::logspace(t_omegas.begin(), t_omegas.end(),
            std::log10(m_omega_range[0]), std::log10(m_omega_range[1]), 10.0_rt);
    }
    else{
        amrex::linspace(t_omegas.begin(), t_omegas.end(),
            m_omega_range[0], m_omega_range[1]);
    }

    m_omegas = amrex::Gpu::DeviceVector<amrex::Real>(m_omega_points);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            t_omegas.begin(), t_omegas.end(), m_omegas.begin());
    amrex::Gpu::Device::streamSynchronize();

    m_has_start     = queryWithParser(pp_radiation, "step_start", m_step_start);
    m_has_stop      = queryWithParser(pp_radiation, "step_stop", m_step_stop);
    m_has_step_skip = queryWithParser(pp_radiation, "step_skip", m_step_skip);
    if (!m_has_step_skip) m_step_skip = 1;

    m_has_start     = queryWithParser(pp_radiation, "step_start", m_step_start);


    if (m_has_start || m_has_stop){
        ablastr::warn_manager::WMRecordWarning(
            "Radiation",
            "Radiation will be integrated from step " + std::to_string(m_step_start) +
            " to step " + std::to_string(m_step_stop),
            ablastr::warn_manager::WarnPriority::low
        );
    }

    if (m_has_step_skip){
        ablastr::warn_manager::WMRecordWarning(
            "Radiation",
            "Radiation.step_skip is set to " + std::to_string(m_step_skip),
            ablastr::warn_manager::WarnPriority::low
        );
    }

    std::vector<std::string> radiation_output_intervals_string = {"0"};
    pp_radiation.queryarr("output_intervals", radiation_output_intervals_string);
    m_output_intervals_parser = utils::parser::IntervalsParser(radiation_output_intervals_string);

    // Cell sizes
    m_d[0] = geom.CellSize(0);
    m_d[1] = geom.CellSize(1);
    m_d[2] = geom.CellSize(2);

    // Shape factor
    m_shape_factor = shape_factor;

    m_center = center;

    m_FF = amrex::Gpu::DeviceVector<amrex::Real>(m_omega_points * det_x.size());

    prepare_data_on_gpu ();

}

void RadiationHandler::prepare_data_on_gpu ()
{
    WARPX_PROFILE("RadiationHandler::prepare_data_on_gpu");

    const auto center = m_center;

    const auto* p_det_x = m_det_x.dataPtr();
    const auto* p_det_y = m_det_y.dataPtr();
    const auto* p_det_z = m_det_z.dataPtr();

    auto* p_det_n_x = m_det_n_x.dataPtr();
    auto* p_det_n_y = m_det_n_y.dataPtr();
    auto* p_det_n_z = m_det_n_z.dataPtr();

    const auto NN = m_det_x.size();
    amrex::ParallelFor(NN, [=] AMREX_GPU_DEVICE (int ip){
        const auto one_over_dist = 1.0_rt/std::sqrt(
            amrex::Math::powi<2>(center[0] - p_det_x[ip]) +
            amrex::Math::powi<2>(center[1] - p_det_y[ip]) +
            amrex::Math::powi<2>(center[2] - p_det_z[ip]));
        p_det_n_x[ip] = (p_det_x[ip]-center[0])*one_over_dist;
        p_det_n_y[ip] = (p_det_y[ip]-center[1])*one_over_dist;
        p_det_n_z[ip] = (p_det_z[ip]-center[2])*one_over_dist;
    });

    const auto dx = m_d[0];
    const auto dy = m_d[1];
    const auto dz = m_d[2];

    constexpr auto inv_c = 1._prt/(PhysConst::c);

    const auto dx_over_2c = dx*inv_c*.5_rt;
    const auto dy_over_2c = dy*inv_c*.5_rt;
    const auto dz_over_2c = dz*inv_c*.5_rt;

    const auto p_omegas = m_omegas.dataPtr();

    auto* p_m_FF = m_FF.dataPtr();

    const auto how_many_det_pos = static_cast<int>(m_det_x.size());
    const auto om_times_det_size =  static_cast<int>(m_omegas.size() * how_many_det_pos);
    amrex::ParallelFor(
        TypeList<CompileTimeOptions<1,2,3>>{},
        {m_shape_factor},
        om_times_det_size, [=] AMREX_GPU_DEVICE (int i_om_det, auto shape_factor_runtime){
            const int i_det = i_om_det % (how_many_det_pos);
            const int i_om  = i_om_det / (how_many_det_pos);

            const auto fcx = p_omegas[i_om]*dx_over_2c;
            const auto fcy = p_omegas[i_om]*dy_over_2c;
            const auto fcz = p_omegas[i_om]*dz_over_2c;

            const auto nx = p_det_n_x[i_det];
            const auto ny = p_det_n_y[i_det];
            const auto nz = p_det_n_z[i_det];

            constexpr auto sinc = [](amrex::Real x){return std::sin(x)/x;};
            const auto FF = amrex::Math::powi<shape_factor_runtime*2>(
                sinc(fcx*nx)*sinc(fcy*ny)*sinc(fcz*nz));

            p_m_FF[i_om_det] = FF;
        });
}


void RadiationHandler::add_radiation_contribution(
    const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc,
    const amrex::Real current_time, const int timestep)
{
    WARPX_PROFILE("RadiationHandler::add_radiation_contribution");

    if (((m_has_start) && (timestep < m_step_start)) ||
        ((m_has_stop) && (timestep > m_step_stop)) ||
        ((m_has_step_skip) && (timestep % m_step_skip != 0))) {
        return;
    }

        constexpr auto c = PhysConst::c;
        constexpr auto inv_c = 1._prt/(PhysConst::c);
        constexpr auto inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);

        for (int lev = 0; lev <= pc->finestLevel(); ++lev)
        {
#ifdef AMREX_USE_OMP
            #pragma omp parallel
#endif
            {
                for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti)
                {

                    long const np = pti.numParticles();
                    const auto& attribs = pti.GetAttribs();
                    const auto* p_w = attribs[PIdx::w].data();
                    const auto* p_ux = attribs[PIdx::ux].data();
                    const auto* p_uy = attribs[PIdx::uy].data();
                    const auto* p_uz = attribs[PIdx::uz].data();

                    const auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                    auto& part = pc->GetParticles(lev)[index];
                    auto& soa = part.GetStructOfArrays();

                    const auto* p_ux_old = soa.GetRealData(pc->GetRealCompIndex("old_u_x")).data();
                    const auto* p_uy_old = soa.GetRealData(pc->GetRealCompIndex("old_u_y")).data();
                    const auto* p_uz_old = soa.GetRealData(pc->GetRealCompIndex("old_u_z")).data();

                    auto GetPosition = GetParticlePosition<PIdx>(pti);
                    auto const q = pc->getCharge();

                    const auto how_many_det_pos = static_cast<int>(m_det_x.size());

                    const auto p_omegas = m_omegas.dataPtr();

                    const auto* p_det_n_x = m_det_n_x.dataPtr();
                    const auto* p_det_n_y = m_det_n_y.dataPtr();
                    const auto* p_det_n_z = m_det_n_z.dataPtr();

                    const auto* p_m_FF = m_FF.dataPtr();

                    auto* p_radiation_data = m_radiation_data.dataPtr();

                    const auto omega_points = m_omega_points;

                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE((np-1) == static_cast<int>(np-1), "too many particles!");


                    m_b_x.resize(np);
                    m_b_y.resize(np);
                    m_b_z.resize(np);
                    m_bp_x.resize(np);
                    m_bp_y.resize(np);
                    m_bp_z.resize(np);

                    auto* p_b_x = m_b_x.dataPtr();
                    auto* p_b_y = m_b_y.dataPtr();
                    auto* p_b_z = m_b_z.dataPtr();
                    auto* p_bp_x = m_bp_x.dataPtr();
                    auto* p_bp_y = m_bp_y.dataPtr();
                    auto* p_bp_z = m_bp_z.dataPtr();


                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip){
                        const auto ux = 0.5_prt*(p_ux[ip] + p_ux_old[ip]);
                        const auto uy = 0.5_prt*(p_uy[ip] + p_uy_old[ip]);
                        const auto uz = 0.5_prt*(p_uz[ip] + p_uz_old[ip]);

                        auto const u2 = ux*ux + uy*uy + uz*uz;

                        auto const one_over_gamma = 1._prt/std::sqrt(1.0_rt + u2*inv_c2);
                        auto const one_over_gamma_c = one_over_gamma*inv_c;
                        const auto one_over_dt_gamma_c = one_over_gamma_c/dt;

                        p_b_x[ip]  = ux*one_over_gamma_c;
                        p_b_y[ip]  = uy*one_over_gamma_c;
                        p_b_z[ip]  = uz*one_over_gamma_c;
                        p_bp_x[ip] = (p_ux[ip] - p_ux_old[ip])*one_over_dt_gamma_c;
                        p_bp_y[ip] = (p_uy[ip] - p_uy_old[ip])*one_over_dt_gamma_c;
                        p_bp_z[ip] = (p_uz[ip] - p_uz_old[ip])*one_over_dt_gamma_c;
                    });

                    const auto np_times_det_pos = np*how_many_det_pos;

                    amrex::ParallelFor(
                        np_times_det_pos, [=] AMREX_GPU_DEVICE(int ii){
                            const int ip  = ii / (np);
                            const int i_det = ii % (np);

                            amrex::ParticleReal xp, yp, zp;
                            GetPosition.AsStored(ip, xp, yp, zp);

                            const auto bx = p_b_x[ip];
                            const auto by = p_b_y[ip];
                            const auto bz = p_b_z[ip];

                            const auto bpx = p_bp_x[ip];
                            const auto bpy = p_bp_y[ip];
                            const auto bpz = p_bp_z[ip];

                            const auto w = p_w[ip];

                            const auto nx = p_det_n_x[i_det];
                            const auto ny = p_det_n_y[i_det];
                            const auto nz = p_det_n_z[i_det];

                            const auto one_minus_b_dot_n = 1.0_prt - (bx*nx + by*ny + bz*nz);

                            //Calculation of 1_beta.n, n corresponds to m_det_direction, the direction of the normal
                            const auto n_minus_beta_x = nx - bx;
                            const auto n_minus_beta_y = ny - by;
                            const auto n_minus_beta_z = nz - bz;

                            //Calculation of nxbeta
                            const auto n_minus_beta_cross_bp_x = n_minus_beta_y*bpz - n_minus_beta_z*bpy;
                            const auto n_minus_beta_cross_bp_y = n_minus_beta_z*bpx - n_minus_beta_x*bpz;
                            const auto n_minus_beta_cross_bp_z = n_minus_beta_x*bpy - n_minus_beta_y*bpx;

                            //Calculation of nxnxbeta
                            const auto n_cross_n_minus_beta_cross_bp_x = ny*n_minus_beta_cross_bp_z - nz*n_minus_beta_cross_bp_y;
                            const auto n_cross_n_minus_beta_cross_bp_y = nz*n_minus_beta_cross_bp_x - nx*n_minus_beta_cross_bp_z;
                            const auto n_cross_n_minus_beta_cross_bp_z = nx*n_minus_beta_cross_bp_y - ny*n_minus_beta_cross_bp_x;

                            const auto nyquist_threshold = ablastr::constant::math::pi/one_minus_b_dot_n/dt;

                            const auto n_dot_r = nx*xp + ny*yp + nz*zp;

#ifdef AMREX_USE_OMP
                            #pragma omp simd
#endif
                            for (int i_om = 0; i_om < omega_points; ++i_om){
                                const auto i_omega_over_c = Complex{0.0_prt, 1.0_prt}*p_omegas[i_om]*inv_c;
                                const auto phase_term = amrex::exp(i_omega_over_c*(c*current_time - (n_dot_r)));

                                const auto FF = p_m_FF[i_om*how_many_det_pos + i_det];
                                const auto form_factor = std::sqrt(w + (w*w-w)*FF);

                                const auto coeff = q*phase_term/(one_minus_b_dot_n*one_minus_b_dot_n)*form_factor;

                                //Nyquist limiter
                                const amrex::Real nyquist_flag = (p_omegas[i_om] < nyquist_threshold);

                                const auto cx = coeff*n_cross_n_minus_beta_cross_bp_x*nyquist_flag;
                                const auto cy = coeff*n_cross_n_minus_beta_cross_bp_y*nyquist_flag;
                                const auto cz = coeff*n_cross_n_minus_beta_cross_bp_z*nyquist_flag;

                                constexpr int ncomp = 3;
                                const int idx0 = (i_om*how_many_det_pos + i_det)*ncomp;
                                const int idx1 = idx0 + 1;
                                const int idx2 = idx0 + 2;

                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx0].m_real, cx.m_real);
                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx0].m_imag, cx.m_imag);
                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx1].m_real, cy.m_real);
                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx1].m_imag, cy.m_imag);
                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx2].m_real, cz.m_real);
                                amrex::HostDevice::Atomic::Add(&p_radiation_data[idx2].m_imag, cz.m_imag);
                            }
                        });
                }
            }
        }
    }

void RadiationHandler::dump_radiation (
    const amrex::Real dt, const int timestep, const std::string& filename)
{
    WARPX_PROFILE("RadiationHandler::dump_radiation");

    if (!m_output_intervals_parser.contains(timestep+1)){ return; }
    Integral_overtime(dt);
    gather_and_write_radiation(filename, timestep);
}

void RadiationHandler::gather_and_write_radiation(const std::string& filename, [[maybe_unused]] const int timestep)
{
    WARPX_PROFILE("RadiationHandler::gather_and_write_radiation");

    auto radiation_data_cpu = amrex::Vector<amrex::Real>(m_det_pts[0]*m_det_pts[1]*m_omega_points);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_radiation_calculation.begin(), m_radiation_calculation.end(), radiation_data_cpu.begin());
    amrex::Gpu::streamSynchronize();

    amrex::ParallelDescriptor::ReduceRealSum(radiation_data_cpu.data(), radiation_data_cpu.size());

    if ( !ParallelDescriptor::IOProcessor() ) { return; }

    const auto how_many = m_det_pts[0]*m_det_pts[1];

    auto det_pos_x_cpu = amrex::Vector<amrex::Real>(how_many);
    auto det_pos_y_cpu = amrex::Vector<amrex::Real>(how_many);
    auto det_pos_z_cpu = amrex::Vector<amrex::Real>(how_many);
    auto omegas_cpu = amrex::Vector<amrex::Real>(m_omega_points);

    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_det_x.begin(), m_det_x.end(), det_pos_x_cpu.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_det_y.begin(), m_det_y.end(), det_pos_y_cpu.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_det_z.begin(), m_det_z.end(), det_pos_z_cpu.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_omegas.begin(), m_omegas.end(), omegas_cpu.begin());
    amrex::Gpu::streamSynchronize();

#ifdef WARPX_USE_OPENPMD
    constexpr int file_min_digits = 6;
    const auto openpmd_backend = WarpXOpenPMDFileType();
    const auto file_suffix = std::string("_%0") + std::to_string(file_min_digits) + std::string("T");
    auto opmd_file_path = filename;
    opmd_file_path.append("/openpmd").append(file_suffix).append(".").append(openpmd_backend);

    // transform paths for Windows
    #ifdef _WIN32
        opmd_file_path = openPMD::auxiliary::replace_all(opmd_file_path, "/", "\\");
    #endif

    // Create the OpenPMD series
    auto series = io::Series(
        opmd_file_path,
        io::Access::CREATE);
    auto i = series.iterations[timestep+1];

    // record
    auto f_rad = i.meshes["radiation_data"];
    auto f_omegas = i.meshes["omegas"];
    auto f_g1 = i.meshes["g1"];
    auto f_g2 = i.meshes["g2"];
    auto f_det_x = i.meshes["det_x"];
    auto f_det_y = i.meshes["det_y"];
    auto f_det_z = i.meshes["det_z"];

    // meta data
    f_rad.setAxisLabels({"i_om", "i_det", "j_det"});
    f_omegas.setAxisLabels({"idx"});
    f_g1.setAxisLabels({"idx"});
    f_g2.setAxisLabels({"idx"});
    f_rad.setGridSpacing<amrex::Real>({1,1,1});
    f_omegas.setGridSpacing<amrex::Real>({1});
    f_g1.setGridSpacing<amrex::Real>({1});
    f_g2.setGridSpacing<amrex::Real>({1});
    f_rad.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt,0.0_rt,0.0_rt});
    f_omegas.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt});
    f_g1.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt});
    f_g2.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt});

    f_det_x.setAxisLabels({"i_det", "j_det"});
    f_det_x.setGridSpacing<amrex::Real>({1,1});
    f_det_x.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt,0.0_rt});
    f_det_y.setAxisLabels({"i_det", "j_det"});
    f_det_y.setGridSpacing<amrex::Real>({1,1});
    f_det_y.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt,0.0_rt});
    f_det_z.setAxisLabels({"i_det", "j_det"});
    f_det_z.setGridSpacing<amrex::Real>({1,1});
    f_det_z.setGridGlobalOffset(std::vector<amrex::Real>{0.0_rt,0.0_rt});

    // record components
    auto rad_data = f_rad[io::RecordComponent::SCALAR];
    auto omegas = f_omegas[io::RecordComponent::SCALAR];
    auto g1s = f_g1[io::RecordComponent::SCALAR];
    auto g2s = f_g2[io::RecordComponent::SCALAR];
    auto det_xs = f_det_x[io::RecordComponent::SCALAR];
    auto det_ys = f_det_y[io::RecordComponent::SCALAR];
    auto det_zs = f_det_z[io::RecordComponent::SCALAR];

    // prepare datasets
    const auto dtype = io::determineDatatype<amrex::Real>();
    rad_data.setPosition<amrex::Real>({0.0,0.0,0.0});
    omegas.setPosition<amrex::Real>({0.0});
    g1s.setPosition<amrex::Real>({0.0});
    g2s.setPosition<amrex::Real>({0.0});
    det_xs.setPosition<amrex::Real>({0.0,0.0});
    det_ys.setPosition<amrex::Real>({0.0,0.0});
    det_zs.setPosition<amrex::Real>({0.0,0.0});
    rad_data.resetDataset(io::Dataset(dtype,
        {(long unsigned int) m_omega_points, (long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]}));
    omegas.resetDataset(io::Dataset(dtype,
        {(long unsigned int) m_omega_points}));
    g1s.resetDataset(io::Dataset(dtype,
        {(long unsigned int) m_det_pts[0]}));
    g2s.resetDataset(io::Dataset(dtype,
        {(long unsigned int) m_det_pts[1]}));
    det_xs.resetDataset(io::Dataset(dtype,
        { (long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]}));
    det_ys.resetDataset(io::Dataset(dtype,
        { (long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]}));
    det_zs.resetDataset(io::Dataset(dtype,
        { (long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]}));

    // write data
    rad_data.storeChunkRaw(
        radiation_data_cpu.data(),
        {0,0,0}, {(long unsigned int) m_omega_points, (long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]});
    omegas.storeChunkRaw(
        omegas_cpu.data(),
        {0}, {(long unsigned int) m_omega_points});
    g1s.storeChunkRaw(
        m_grid[0].data(),
        {0}, {(long unsigned int) m_grid[0].size()});
    g2s.storeChunkRaw(
        m_grid[1].data(),
        {0}, {(long unsigned int) m_grid[1].size()});
    det_xs.storeChunkRaw(
        det_pos_x_cpu.data(),
        {0,0}, {(long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]});
    det_ys.storeChunkRaw(
        det_pos_y_cpu.data(),
        {0,0}, {(long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]});
    det_zs.storeChunkRaw(
        det_pos_z_cpu.data(),
        {0,0}, {(long unsigned int) m_det_pts[0], (long unsigned int) m_det_pts[1]});
    series.flush();
#else

    auto of = std::ofstream(filename, std::ios::binary);

    int idx = 0;
    for(int i_om=0; i_om < m_omega_points; ++i_om){
        for (int i_det = 0; i_det < how_many; ++i_det)
        {
#if defined(WARPX_DIM_3D)
            of << omegas_cpu[i_om] << " "  << det_pos_x_cpu[i_det] << " " << det_pos_y_cpu[i_det] << " " << det_pos_z_cpu[i_det]  << " " << radiation_data_cpu[idx++] << "\n";
#elif defined(WARPX_DIM_XZ)
            of << omegas_cpu[i_om] << " " << det_pos_x_cpu[i_det] << " " << det_pos_z_cpu[i_det]  << " " << radiation_data_cpu[idx++] << "\n";
#endif
        }
    }

    of.close();

#endif

}

void RadiationHandler::Integral_overtime(const amrex::Real dt)
{
    WARPX_PROFILE("RadiationHandler::Integral_overtime");

    const amrex::Real long_dt = dt * m_step_skip;

    const auto factor = long_dt*long_dt/16/std::pow(ablastr::constant::math::pi,3)/PhysConst::ep0/(PhysConst::c);

    const auto how_many = m_det_pts[0]*m_det_pts[1];

    auto* const p_radiation_data = m_radiation_data.dataPtr();

    m_radiation_calculation = amrex::Gpu::DeviceVector<amrex::Real>(how_many*m_omega_points);
    auto* const p_radiation_calculation = m_radiation_calculation.dataPtr();

    amrex::ParallelFor(m_omega_points*how_many,
        [=] AMREX_GPU_DEVICE(int idx){
            const int idx0 = idx*3;
            const int idx1 = idx0 + 1;
            const int idx2 = idx0 + 2;
            p_radiation_calculation[idx]=(amrex::norm(p_radiation_data[idx0]) + amrex::norm(p_radiation_data[idx1]) + amrex::norm(p_radiation_data[idx2]))*factor;

        });
}
