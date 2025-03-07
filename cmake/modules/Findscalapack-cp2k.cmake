if(NOT cp2k::ScaLAPACK)
  add_library(cp2k::ScaLAPACK INTERFACE IMPORTED)
endif()

if(CP2K_SCALAPACK_VENDOR STREQUAL "Intel" OR CP2K_SCALAPACK_VENDOR STREQUAL "MKL")
  # Use MKL config file to find ScaLAPACK
  find_package(onemkl-cp2k)
  target_link_libraries(cp2k::ScaLAPACK INTERFACE MKL::MKL_BLACS MKL::MKL_SCALAPACK)
else()
  # Use CP2K FindSCALAPACK.cmake to find ScaLAPACK
  find_package(SCALAPACK)
endif()
