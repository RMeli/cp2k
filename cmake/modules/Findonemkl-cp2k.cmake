# TODO: Define at CP2K level?
if(CP2K_USE_BLAS_STATIC)
  set(MKL_STATIC "static")
else()
  set(MKL_STATIC "dynamic")
endif()

if(CP2K_BLAS_INTERFACE STREQUAL "32bits")
  set(MKL_INTERFACE "lp64")
else()
  set(MKL_INTERFACE "ilp64")
endif()

if(CP2K_BLAS_THREADING STREQUAL "gnu-thread")
  set(MKL_THREADING "gnu_thread")
elseif(CP2K_BLAS_THREADING STREQUAL "intel-thread")
  set(MKL_THREADING "intel_thread")
elseif(CP2K_BLAS_THREADING STREQUAL "tbb-thread")
  set(MKL_THREADING "tbb_thread")
elseif(CP2K_BLAS_THREADING STREQUAL "sequential")
  set(MKL_THREADING "sequential")
endif()

# TODO: Cover OpenMPI
set(MKL_MPI "mpich")

find_package(MKL REQUIRED CONFIG)

get_target_property(CP2K_BLAS_LIBRARIES MKL::MKL INTERFACE_LINK_LIBRARIES)
get_target_property(CP2K_LAPACK_LIBRARIES MKL::MKL INTERFACE_LINK_LIBRARIES)

# TODO: Support old versions of MKL without CMake config?
