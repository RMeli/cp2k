# Try to find NVPL with the provided CONFIG
find_package(nvpl CONFIG COMPONENTS blac lapack scalapack fft)

if(CP2K_BLAS_INTERFACE STREQUAL "64bits")
  set(_nvpl_int "_ilp64")
else()
  set(_nvpl_int "_lp64")
endif()

if(CP2K_BLAS_THREADING STREQUAL "sequential")
  set(_nvpl_thread "_seq")
else()
  set(_nvpl_thread "_omp")
endif()

# Look for separate components
if(NOT nvpl_FOUND),
  find_package(nvpl_blas CONFIG REQUIRED)
  find_package(nvpl_lapack CONFIG REQUIRED)
  find_package(nvpl_scalapack CONFIG REQUIRED)
  # 0.4.0 is the first version with the fftw3.h header
  find_package(nvpl_fftw 0.4.0 CONFIG REQUIRED)
endif()

if(NOT TARGET "cp2k::BLAS::NVPL::blas")
  add_library("cp2k::BLAS::NVPL::blas" INTERFACE IMPORTED)
  target_link_libraries("cp2k::BLAS::NVPL::blas" INTERFACE "nvpl::blas${_nvpl_int}${_nvpl_thread}")
endif()

if(NOT TARGET "cp2k::BLAS::NVPL::lapack")
  add_library("cp2k::BLAS::NVPL::lapack" INTERFACE IMPORTED)
  target_link_libraries("cp2k::BLAS::NVPL::lapack" INTERFACE "nvpl::lapack${_nvpl_int}${_nvpl_thread}")
endif()

if(NOT TARGET "cp2k::BLAS::NVPL::scalapack_link")
  add_library("cp2k::BLAS::NVPL::scalapack_link" INTERFACE IMPORTED)

  get_target_property(CP2K_NVPL_LAPACK_LIBRARIES "nvpl::lapack${_nvpl_int}${_nvpl_thread}" INTERFACE_LINK_LIBRARIES)
  get_target_property(CP2K_NVPL_SCALAPACK_LIBRARIES "nvpl::scalapack${_nvpl_int}" INTERFACE_LINK_LIBRARIES)
  get_target_property(CP2K_NVPL_BLAS_INCLUDE_DIRS "nvpl::blas${_nvpl_int}${_nvpl_thread}" INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_NVPL_LAPACK_INCLUDE_DIRS "nvpl::lapack${_nvpl_int}${_nvpl_thread}" INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_NVPL_SCALAPACK_INCLUDE_DIRS "nvpl::scalapack${_nvpl_int}" INTERFACE_INCLUDE_DIRECTORIES)

  set_target_properties(
    cosma::BLAS::NVPL::scalapack_link 
    PROPERTIES INTERFACE_LINK_LIBRARIES 
    "${COSMA_NVPL_LAPACK_LIBRARIES};${COSMA_NVPL_SCALAPACK_LIBRARIES}")
  set_target_properties(
    cosma::BLAS::NVPL::scalapack_link 
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES 
    "${COSMA_NVPL_BLAS_INCLUDE_DIRS};${COSMA_NVPL_LAPACK_INCLUDE_DIRS};${COSMA_NVPL_SCALAPACK_INCLUDE_DIRS}")
endif()

if(NOT TARGET "cp2k::FFTW3::fftw3")
  add_library("cp2k::FFTW3::fftw3" INTERFACE IMPORTED)
  target_link_libraries("cp2k::FFT3W::fftw3" INTERFACE "nvpl::fftw")
endif()
