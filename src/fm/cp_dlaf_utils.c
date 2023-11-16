/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2023 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__DLAF)

#include <dlaf_c/eigensolver/eigensolver.h>
#include <dlaf_c/factorization/cholesky.h>
#include <dlaf_c/init.h>

void dlaf_init() {
  const char *pika_argv[] = {"cp2k", "--pika:print-bind"};
  const char *dlaf_argv[] = {"cp2k"};
  dlaf_initialize(2, pika_argv, 1, dlaf_argv);
}

void dlaf_pdsyevd_wrapper(int n, double *a, int ia, int ja, int desca[9],
                          double *w, double *z, int iz, int jz, int descz[9],
                          int *info) {
  char uplo = 'L'; // Only lower triangular matrices are supported
  dlaf_pdsyevd(uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, info);
}

void dlaf_pdpotrf_wrapper(int uplo, int n, double *a, int ia, int ja,
                          int desca[9], int *info) {
  char _uplo = ' ';
  if (uplo == 0) {
    _uplo = 'L';
  } else if (uplo == 1) {
    _uplo = 'U';
  }
  dlaf_pdpotrf(_uplo, n, a, ia, ja, desca, info);
}

void dlaf_pspotrf_wrapper(int uplo, int n, float *a, int ia, int ja,
                          int desca[9], int *info) {
  char _uplo = ' ';
  if (uplo == 0) {
    _uplo = 'L';
  } else if (uplo == 1) {
    _uplo = 'U';
  }
  dlaf_pspotrf(_uplo, n, a, ia, ja, desca, info);
}

#endif // __DLAF
