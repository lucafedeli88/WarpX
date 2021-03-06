/* Copyright 2019-2020 Axel Huebl, Junmin Gu, Maxence Thevenet
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_OPEN_PMD_H_
#define WARPX_OPEN_PMD_H_

#include "Particles/WarpXParticleContainer.H"
#include "Particles/MultiParticleContainer.H" // PIdx
#include "Particles/PhysicalParticleContainer.H"

#include "Diagnostics/ParticleDiag/ParticleDiag.H"

#include <AMReX_AmrParticles.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

#include <memory>
#include <string>
#include <vector>


//
// helper class
//
class Timer
{
public:
  Timer(const char* tag) {m_Tag = tag; m_Start = amrex::second();}
  ~Timer() {
      m_End = amrex::second();
      amrex::ParallelDescriptor::ReduceRealMax(m_End, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print()<<m_Tag<<" took:"<<m_End - m_Start<<" seconds\n";
  }
private:
  amrex::Real m_Start;
  amrex::Real m_End;
  std::string m_Tag;
};


//
//
class WarpXParticleCounter
{
public:
  using ParticleContainer = typename amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;
  using ParticleIter = typename amrex::ParIter<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;

  WarpXParticleCounter(ParticleContainer* pc);
  unsigned long GetTotalNumParticles() {return m_Total;}

  std::vector<unsigned long long> m_ParticleOffsetAtRank;
  std::vector<unsigned long long> m_ParticleSizeAtRank;
private:
  /** get the offset in the overall particle id collection
  *
  * @param[out] numParticles particles on this processor  / amrex fab
  * @param[out] offset particle offset over all, mpi-global amrex fabs
  * @param[out] sum number of all particles from all amrex fabs
  */
  void GetParticleOffsetOfProcessor(const long& numParticles,
                    unsigned long long& offset,
                    unsigned long long& sum)  const ;


  int m_MPIRank = 0;
  int m_MPISize = 1;

  unsigned long long m_Total = 0;

  std::vector<unsigned long long> m_ParticleCounterByLevel;
};


#ifdef WARPX_USE_OPENPMD
//
//
/** Writer logic for openPMD particles and fields */
class WarpXOpenPMDPlot
{
public:
  using ParticleContainer = typename amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;
  using ParticleIter = typename amrex::ParConstIter<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;

  /** Initialize openPMD I/O routines
   *
   * @param oneFilePerTS write one file per timestep
   * @param filetype file backend, e.g. "bp" or "h5"
   * @param fieldPMLdirections PML field solver, @see WarpX::getPMLdirections()
   */
  WarpXOpenPMDPlot (bool oneFilePerTS, std::string filetype, std::vector<bool> fieldPMLdirections);

  ~WarpXOpenPMDPlot ();

  /** Set Iteration Step for the series
   *
   * @note If an iteration has been written, then it will give a warning
   *
   */
  void SetStep (int ts, const std::string& dirPrefix,
                bool isBTD=false);

  /** Close the step
   *
   * Signal that no further updates will be written for the step.
   */
  void CloseStep (bool isBTD = false, bool isLastBTDFlush = false);

  void WriteOpenPMDParticles (const amrex::Vector<ParticleDiag>& particle_diags);

  void WriteOpenPMDFields (
              const std::vector<std::string>& varnames,
              const amrex::MultiFab& mf,
              const amrex::Geometry& geom,
              const int iteration, const double time,
              bool isBTD = false,
              const amrex::Geometry& full_BTD_snapshot=amrex::Geometry() ) const;


private:
  void Init (openPMD::Access access, bool isBTD);

  /** This function sets up the entries for storing the particle positions, global IDs, and constant records (charge, mass)
  *
  * @param[in] currSpecies Corresponding openPMD species
  * @param[in] np          Number of particles
  * @param[in] charge      Charge of the particles (note: fix for ions)
  * @param[in] mass        Mass of the particles
  */
  void SetupPos (
        openPMD::ParticleSpecies& currSpecies,
        const unsigned long long& np,
        amrex::ParticleReal const charge,
        amrex::ParticleReal const mass) const ;

  /** This function sets up the entries for particle properties
   *
   * @param[in] currSpecies The openPMD species
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   * @param[in] write_int_comp The int attribute ids, from WarpX
   * @param[in] int_comp_names The int attribute names, from WarpX
   * @param[in] np  Number of particles
   */
  void SetupRealProperties (openPMD::ParticleSpecies& currSpecies,
               const amrex::Vector<int>& write_real_comp,
               const amrex::Vector<std::string>& real_comp_names,
               const amrex::Vector<int>& write_int_comp,
               const amrex::Vector<std::string>& int_comp_names,
               unsigned long long np) const;

  /** This function saves the values of the entries for particle properties
   *
   * @param[in] pti WarpX particle iterator
   * @param[in] currSpecies The openPMD species to save to
   * @param[in] offset offset to start saving  the particle iterator contents
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   * @param[in] write_int_comp The int attribute ids, from WarpX
   * @param[in] int_comp_names The int attribute names, from WarpX
   */
  void SaveRealProperty (ParticleIter& pti, //int, int,
            openPMD::ParticleSpecies& currSpecies,
            unsigned long long offset,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<std::string>& real_comp_names,
            const amrex::Vector<int>& write_int_comp,
            const amrex::Vector<std::string>& int_comp_names) const;

  /** This function saves the plot file
   *
   * @param[in] pc WarpX particle container
   * @param[in] name species name
   * @param[in] iteration timestep
   * @param[in] write_real_comp The real attribute ids, from WarpX
   * @param[in] real_comp_names The real attribute names, from WarpX
   * @param[in] write_int_comp The int attribute ids, from WarpX
   * @param[in] int_comp_names The int attribute names, from WarpX
   * @param[in] charge         Charge of the particles (note: fix for ions)
   * @param[in] mass           Mass of the particles
   */
  void DumpToFile (ParticleContainer* pc,
            const std::string& name,
            int iteration,
            const amrex::Vector<int>& write_real_comp,
            const amrex::Vector<int>& write_int_comp,
            const amrex::Vector<std::string>& real_comp_names,
            const amrex::Vector<std::string>&  int_comp_names,
            amrex::ParticleReal const charge,
            amrex::ParticleReal const mass) const;

  /** Get the openPMD-api filename for openPMD::Series
   *
   * No need for ts in the file name, openPMD handles steps (iterations).
   *
   * @param[inout] filepath the path and filename for openPMD::Series
   *               passes a prefix path in and appends the filename
   * @return pure filename w/o path
   */
  std::string GetFileName (std::string& filepath);

  std::unique_ptr<openPMD::Series> m_Series;

  /** This is the output directory
   *
   * This usually does not yet end in a `/`.
   * It does not yet include the file prefix of the openPMD series, which will
   * be appended by the GetFileName function.
   */
  std::string m_dirPrefix;

  int m_MPIRank = 0;
  int m_MPISize = 1;

  int m_NumSoARealAttributes = PIdx::nattribs; //! WarpX' additional real particle attributes in SoA
  int m_NumAoSRealAttributes = 0; //! WarpX definition: no additional real attributes in particle AoS

  //int m_NumSoAIntAttributes = PIdx::nattribs; //! WarpX' additional int particle attributes in SoA
  int m_NumAoSIntAttributes = 0; //! WarpX definition: no additional int attributes in particle AoS

  bool m_OneFilePerTS = true;  //! write in openPMD fileBased manner for individual time steps
  std::string m_OpenPMDFileType = "bp"; //! MPI-parallel openPMD backend: bp or h5
  int m_CurrentStep  = -1;

  // meta data
  std::vector< bool > m_fieldPMLdirections; //! @see WarpX::getPMLdirections()
};
#endif // WARPX_USE_OPENPMD

#endif // WARPX_OPEN_PMD_H
