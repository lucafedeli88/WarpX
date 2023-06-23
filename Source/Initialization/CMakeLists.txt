foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(WarpX_${SD}
      PRIVATE
        WarpXInit.cpp
        InjectorDensity.cpp
        InjectorMomentum.cpp
        PlasmaInjector.cpp
        WarpXInitData.cpp
        TemperatureProperties.cpp
        VelocityProperties.cpp
        GetTemperature.cpp
        GetVelocity.cpp
    )
endforeach()
