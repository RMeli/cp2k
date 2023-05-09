//
// Distributed Linear Algebra with Future (DLAF)
//
// Copyright (c) 2018-2023, ETH Zurich
// All rights reserved.
//
// Please, refer to the LICENSE file in the root directory.
// SPDX-License-Identifier: BSD-3-Clause
//

#include <blas/util.hh>
#include <cstdlib>

#include <dlaf/communication/communicator.h>
#include <dlaf/communication/communicator_grid.h>
#include <dlaf/communication/error.h>
#include <dlaf/eigensolver/eigensolver.h>
#include <dlaf/factorization/cholesky.h>
#include <dlaf/init.h>
#include <dlaf/matrix/distribution.h>
#include <dlaf/matrix/index.h>
#include <dlaf/matrix/layout_info.h>
#include <dlaf/matrix/matrix.h>
#include <dlaf/matrix/matrix_mirror.h>
#include <dlaf/types.h>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <pika/execution.hpp>
#include <pika/init.hpp>
#include <pika/program_options.hpp>
#include <pika/runtime.hpp>

// TODO: Remove once https://github.com/eth-cscs/DLA-Future/pull/668 is
// merged.
// #include <mkl_service.h>
// #include <omp.h>

static bool dlaf_init_ = false;

// Cblacs does not seem to have public headers files and is only used in
// scalapack.

extern "C" MPI_Comm Cblacs2sys_handle(int ictxt);
extern "C" void Cblacs_get(int ictxt, int inum, int *comm);
extern "C" void Cblacs_gridinfo(int ictxt, int *np, int *mp, int *px, int *py);

// queries the grid blacs context to get the communication blacs context
static int get_comm_context(const int grid_context) {
  int comm_context;
  Cblacs_get(grid_context, 10, &comm_context);
  return comm_context;
}

// gets MPI_Comm from the grid blacs context
static MPI_Comm get_communicator(const int grid_context) {
  int comm_context = get_comm_context(grid_context);
  MPI_Comm comm = Cblacs2sys_handle(comm_context);
  return comm;
}

static int get_grid_context(int *desca) { return desca[1]; }

extern "C" void dlaf_init() {
  if (!dlaf_init_) {
    int argc = 1;
    const char *const argv[] = {"cp2k", nullptr};

    pika::program_options::options_description desc("cp2k");
    desc.add(dlaf::getOptionsDescription());

    /* pika initialization */
    pika::init_params p;
    p.rp_callback = dlaf::initResourcePartitionerHandler;
    p.desc_cmdline = desc;
    pika::start(nullptr, argc, argv, p);

    /* DLA-Future initialization */
    dlaf::initialize(argc, argv);
    dlaf_init_ = true;
    pika::suspend();
  }
}

extern "C" void dlaf_finalize() {
  pika::resume();
  pika::finalize();
  dlaf::finalize();
  pika::stop();
  dlaf_init_ = false;
}

// TODO: Remove
// class single_threaded_omp {
// public:
//   single_threaded_omp() : old_threads(mkl_get_max_threads()) {
//     mkl_set_num_threads(1);
//   }
//   ~single_threaded_omp() { mkl_set_num_threads(old_threads); }
//
// private:
//   int old_threads:;
// };

std::tuple<dlaf::matrix::Distribution, dlaf::matrix::LayoutInfo,
           dlaf::comm::CommunicatorGrid>
dlaf_setup(int *desca__) {
  int m, n;   // Matrix sizes
  int nb, mb; // Block sizes

  // retrive the matrix sizes
  m = desca__[3];
  n = desca__[2];

  nb = desca__[5];
  mb = desca__[4];

  MPI_Comm comm = get_communicator(desca__[1]);

  dlaf::comm::Communicator world(comm);
  DLAF_MPI_CHECK_ERROR(MPI_Barrier(world));

  int dims[2] = {0, 0};
  int periods[2] = {0, 0};
  int coords[2] = {-1, -1};

  Cblacs_gridinfo(desca__[1], dims, dims + 1, coords, coords + 1);

  dlaf::comm::CommunicatorGrid comm_grid(world, dims[0], dims[1],
                                         dlaf::common::Ordering::RowMajor);

  dlaf::GlobalElementSize matrix_size(n, m);
  dlaf::TileElementSize block_size(nb, mb);
  dlaf::comm::Index2D src_rank_index(0, 0); // TODO: Get from BLACS?

  dlaf::matrix::Distribution distribution(matrix_size, block_size,
                                          comm_grid.size(), comm_grid.rank(),
                                          src_rank_index);

  const int ld_ = desca__[8]; // Leading dimesnion
  dlaf::matrix::LayoutInfo layout = colMajorLayout(distribution, ld_);

  return std::make_tuple(distribution, layout, comm_grid);
}

void dlaf_check(char uplo, int *desca, int &info) {
  // TODO: Figure out why this segfaults for the Eigensolver
  if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l') {
    std::cerr << "ERROR: The UpLo parameter has a incorrect value: '" << uplo;
    std::cerr << "'. Please check the ScaLAPACK documentation.\n";
    info = -1;
    return;
  }

  if (desca[0] != 1) {
    info = -1;
    std::cerr << "ERROR: DLA-Future can only treat dense matrices.\n";
    return;
  }

  if (!dlaf_init_) {
    std::cerr << "Error: DLA-Future must be initialized.\n";
    info = -1;
    return;
  }
}

template <typename T>
void pxpotrf_dlaf(char uplo__, int n__, T *a__, int ia__, int ja__,
                  int *desca__, int &info__) {

  dlaf_check(uplo__, desca__, info__);
  if (info__ == -1)
    return;

  // TODO: Remove
  //  single_threaded_omp sto{};

  pika::resume();

  // TODO:
  // DONE - dlaf initialization
  // DONE - matrix mirror
  // DONE - resume suspend pika runtime
  //      - general cleanup
  // DONE - fortran interface uplo
  //      - cblcs call
  // DONE - remove omp/mkl_set_num_threads calls (will be handled by DLAF)

  auto [distribution, layout, comm_grid] = dlaf_setup(desca__);

  dlaf::matrix::Matrix<T, dlaf::Device::CPU> mat(std::move(distribution),
                                                 layout, a__);

  {
    dlaf::matrix::MatrixMirror<T, dlaf::Device::Default, dlaf::Device::CPU>
        matrix(mat);

    // uplo__ checked above
    auto dlaf_uplo =
        uplo__ == 'U' or uplo__ == 'u' ? blas::Uplo::Upper : blas::Uplo::Lower;

    dlaf::factorization::cholesky<dlaf::Backend::Default, dlaf::Device::Default,
                                  T>(comm_grid, dlaf_uplo, matrix.get());
  } // Destroy MatrixMirror; copy results back to Device::CPU

  mat.waitLocalTiles();

  pika::suspend();
  info__ = 0;
}

extern "C" void pdpotrf_dlaf(char *uplo__, int n__, double *a__, int ia__,
                             int ja__, int *desca__, int *info__) {
  pxpotrf_dlaf<double>(*uplo__, n__, a__, ia__, ja__, desca__, *info__);
}

extern "C" void pspotrf_dlaf(char *uplo__, int n__, float *a__, int ia__,
                             int ja__, int *desca__, int *info__) {
  pxpotrf_dlaf<float>(*uplo__, n__, a__, ia__, ja__, desca__, *info__);
}

template <typename T>
void pxsyevd_dlaf(char jobz__, char uplo__, int n__, T *a__, int ia__, int ja__,
                  int *desca__, T *w__, T *z__, int iz__, int jz__,
                  int *desc_z__, T *work__, int lwork__, int *iwork__,
                  int liwork__, int &info__) {
  dlaf_check(uplo__, desca__, info__);
  if (info__ == -1)
    return;

  // TODO: Remove
  //  single_threaded_omp sto{};

  pika::resume();

  // TODO: Remove
  int rank = 0;
  MPI_Comm comm = get_communicator(desca__[1]);
  MPI_Comm_rank(comm, &rank);

  auto [distribution, layout, comm_grid] = dlaf_setup(desca__);

  // DLAF matrix, manages the correct dependencies of taks involving tiles
  dlaf::matrix::Matrix<T, dlaf::Device::CPU> host_matrix(distribution, layout,
                                                         a__);

  // uplo__ checked above
  auto dlaf_uplo =
      uplo__ == 'U' or uplo__ == 'u' ? blas::Uplo::Upper : blas::Uplo::Lower;

  // TODO: Remove
  if (rank == 0)
    std::cerr << "Calling DLAF eigensolver..." << std::endl;

  // Create DLAF matrices from CP2K-allocated ones
  dlaf::matrix::Matrix<T, dlaf::Device::CPU> eigenvectors_cp2k(
      distribution, layout, z__); // Distributed eigenvectors
  auto eigenvalues_cp2k =
      dlaf::matrix::createMatrixFromColMajor<dlaf::Device::CPU>(
          {n__, 1}, {distribution.blockSize().rows(), 1}, n__, w__);

  {
    // Create matrix mirrors
    dlaf::matrix::MatrixMirror<T, dlaf::Device::Default, dlaf::Device::CPU>
        matrix(host_matrix);
    dlaf::matrix::MatrixMirror<T, dlaf::Device::Default, dlaf::Device::CPU>
        eigenvalues(eigenvalues_cp2k);
    dlaf::matrix::MatrixMirror<T, dlaf::Device::Default, dlaf::Device::CPU>
        eigenvectors(eigenvectors_cp2k);

    // WARN: Hard-coded to LOWER, use dlaf_uplo instead
    dlaf::eigensolver::eigensolver<dlaf::Backend::Default,
                                   dlaf::Device::Default, T>(
        comm_grid, blas::Uplo::Lower, matrix.get(), eigenvalues.get(),
        eigenvectors.get());
  } // Destroy mirrors

  eigenvalues_cp2k.waitLocalTiles();

  // TODO: Remove
  if (rank == 0)
    std::cerr << "DLAF eigensolver terminated successfully!" << std::endl;

  pika::suspend();
  info__ = 0;
}

extern "C" void pdsyevd_dlaf(char *jobz__, char *uplo__, int n__, double *a__,
                             int ia__, int ja__, int *desca__, double *w__,
                             double *z__, int iz__, int jz__, int *desc_z__,
                             double *work__, int lwork__, int *iwork__,
                             int liwork__, int *info__) {
  pxsyevd_dlaf<double>(*jobz__, *uplo__, n__, a__, ia__, ja__, desca__, w__,
                       z__, iz__, jz__, desc_z__, work__, lwork__, iwork__,
                       liwork__, *info__);
}

extern "C" void pssyevd_dlaf(char *jobz__, char *uplo__, int n__, float *a__,
                             int ia__, int ja__, int *desca__, float *w__,
                             float *z__, int iz__, int jz__, int *desc_z__,
                             float *work__, int lwork__, int *iwork__,
                             int liwork__, int *info__) {
  pxsyevd_dlaf<float>(*jobz__, *uplo__, n__, a__, ia__, ja__, desca__, w__, z__,
                      iz__, jz__, desc_z__, work__, lwork__, iwork__, liwork__,
                      *info__);
}
