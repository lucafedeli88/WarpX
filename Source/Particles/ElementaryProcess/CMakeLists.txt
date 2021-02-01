target_sources(WarpX
  PRIVATE
    Ionization.cpp
)

if(WarpX_QED)
    target_sources(WarpX
      PRIVATE
        QEDPairGeneration.cpp
        QEDPhotonEmission.cpp
        QedSpeciesModule.cpp
    )
    add_subdirectory(QEDInternals)
endif()
