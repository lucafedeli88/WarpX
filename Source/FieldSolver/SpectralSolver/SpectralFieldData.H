/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_FIELD_DATA_H_
#define WARPX_SPECTRAL_FIELD_DATA_H_

#include "Utils/WarpX_Complex.H"
#include "SpectralKSpace.H"
#include "AnyFFT.H"

#include <AMReX_MultiFab.H>

#include <string>

// Declare type for spectral fields
using SpectralField = amrex::FabArray< amrex::BaseFab <Complex> >;

/** Index for the regular fields, when stored in spectral space:
 *  - n_fields is automatically the total number of fields
 *  - divE reuses the memory slot for Bx, since Bx is not used when computing divE
 */
struct SpectralFieldIndex {
  enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new, n_fields, divE=3 };
};

/* Index for the regular fields + averaged fields, when stored in spectral space */
struct SpectralAvgFieldIndex {
  enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new, Ex_avg, Ey_avg, Ez_avg, Bx_avg, By_avg, Bz_avg,n_fields };
  // n_fields is automatically the total number of fields
};

/** Index for the PML fields, when stored in spectral space,
 *  (n_fields is automatically set to the total number of fields)
 *  TODO How to include the diagonal components only when needed?
 */
struct SpectralPMLIndex {
  enum {Exx=0, Exy, Exz, Eyx, Eyy, Eyz, Ezx, Ezy, Ezz,
        Bxx  , Bxy, Bxz, Byx, Byy, Byz, Bzx, Bzy, Bzz, Fx, Fy, Fz, Gx, Gy, Gz, n_fields};
};

/** \brief Class that stores the fields in spectral space, and performs the
 *  Fourier transforms between real space and spectral space
 */
class SpectralFieldData
{

    public:
        SpectralFieldData( const int lev,
                           const amrex::BoxArray& realspace_ba,
                           const SpectralKSpace& k_space,
                           const amrex::DistributionMapping& dm,
                           const int n_field_required,
                           const bool periodic_single_box);
        SpectralFieldData() = default; // Default constructor
        SpectralFieldData& operator=(SpectralFieldData&& field_data) = default;
        ~SpectralFieldData();

        void ForwardTransform (const int lev,
                               const amrex::MultiFab& mf, const int field_index,
                               const int i_comp, const amrex::IntVect& stag);
        AMREX_FORCE_INLINE
        void ForwardTransform (const int lev,
                               const amrex::MultiFab& mf, const int field_index, const int i_comp)
        {
            ForwardTransform(lev, mf, field_index, i_comp, mf.ixType().toIntVect());
        }

        void BackwardTransform (const int lev, amrex::MultiFab& mf, const int field_index, const int i_comp);

        // `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

    private:
        // tmpRealField and tmpSpectralField store fields
        // right before/after the Fourier transform
        SpectralField tmpSpectralField; // contains Complexs
        amrex::MultiFab tmpRealField; // contains Reals
        AnyFFT::FFTplans forward_plan, backward_plan;
        // Correcting "shift" factors when performing FFT from/to
        // a cell-centered grid in real space, instead of a nodal grid
        SpectralShiftFactor xshift_FFTfromCell, xshift_FFTtoCell,
                            zshift_FFTfromCell, zshift_FFTtoCell;
#if (AMREX_SPACEDIM==3)
        SpectralShiftFactor yshift_FFTfromCell, yshift_FFTtoCell;
#endif

        bool m_periodic_single_box;
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
