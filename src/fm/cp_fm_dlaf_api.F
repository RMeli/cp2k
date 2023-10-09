!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_fm_dlaf_api

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_CHAR, C_FLOAT, C_LOC

   USE cp_fm_diag_utils, ONLY: cp_fm_redistribute_start, cp_fm_redistribute_end
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_types, ONLY: cp_fm_type
   USE cp_fm_basic_linalg, ONLY: cp_fm_upper_to_full
   USE kinds, ONLY: dp

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_dlaf_api'

#if defined(__DLAF)
   PUBLIC :: cp_pdpotrf_dlaf, cp_pspotrf_dlaf
#endif
   PUBLIC :: cp_fm_diag_dlaf

CONTAINS

#if defined(__DLAF)

!***************************************************************************************************
!> \brief Cholesky factorization using DLA-Future
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
!***************************************************************************************************
   SUBROUTINE cp_pdpotrf_dlaf(uplo, n, a, ia, ja, desca, info)
      CHARACTER, INTENT(IN) ::           uplo
      INTEGER, INTENT(IN)       ::       ia, ja, n
      INTEGER, TARGET     ::      info
      INTEGER, dimension(9)       ::        DESCA
      REAL(KIND=dp), DIMENSION(:, :), TARGET   :: A

      INTERFACE
         SUBROUTINE pdpotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pdpotrf_wrapper')
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            INTEGER(KIND=C_INT), VALUE :: UPLO
            INTEGER(KIND=C_INT), VALUE      ::        IA, JA, N
            TYPE(C_PTR), VALUE ::       INFO
            INTEGER(kind=C_INT), DIMENSION(*)         :: DESCA
            TYPE(C_PTR), VALUE ::              A
         END SUBROUTINE pdpotrf_dlaf
      END INTERFACE

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(A))
#endif

      IF (uplo == "L") THEN
         CALL pdpotrf_dlaf(0, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      ELSE IF (uplo == "U") THEN
         CALL pdpotrf_dlaf(1, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      END IF
   END SUBROUTINE cp_pdpotrf_dlaf

!***************************************************************************************************
!> \brief Cholesky factorization using DLA-Future
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
!***************************************************************************************************
   SUBROUTINE cp_pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO)
      CHARACTER, INTENT(IN) ::           UPLO
      INTEGER, INTENT(IN)       ::       IA, JA, N
      INTEGER, TARGET     ::      INFO
      INTEGER, dimension(9)      ::        DESCA
      REAL, dimension(:, :), TARGET    ::              A

      INTERFACE
         SUBROUTINE pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pspotrf_wrapper')
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            INTEGER(KIND=C_INT), VALUE :: UPLO
            INTEGER(KIND=C_INT), VALUE      ::        IA, JA, N
            TYPE(C_PTR), VALUE    ::       INFO
            INTEGER(kind=C_INT), DIMENSION(*)         :: DESCA
            TYPE(C_PTR), VALUE   ::              A
         END SUBROUTINE pspotrf_dlaf
      END INTERFACE

      IF (uplo == "L") THEN
         CALL pspotrf_dlaf(0, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      ELSE IF (uplo == "U") THEN
         CALL pspotrf_dlaf(1, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      END IF
   END SUBROUTINE cp_pspotrf_dlaf

#endif

   SUBROUTINE cp_fm_diag_dlaf(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), INTENT(IN)  :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues

      CHARACTER(LEN=*), PARAMETER              :: routineN = 'cp_fm_syevd'

      INTEGER                                  :: handle, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET :: eig
      TYPE(cp_fm_type)                         :: eigenvectors_new, matrix_new
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m

      CALL timeset(routineN, handle)

      n = matrix%matrix_struct%nrow_global
      ALLOCATE (eig(n), source=-42.0_dp)

      CALL cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new)

      ! Call scalapack on CPUs that hold the new matrix
      IF (ASSOCIATED(matrix_new%matrix_struct)) THEN
         CALL cp_fm_diag_dlaf_base(matrix_new, eigenvectors_new, eig)
      END IF

      CALL cp_fm_redistribute_end(matrix, eigenvectors, eig, matrix_new, eigenvectors_new)

      nmo = SIZE(eigenvalues, 1)
      IF (nmo > n) THEN
         eigenvalues(1:n) = eig(1:n)
      ELSE
         eigenvalues(1:nmo) = eig(1:nmo)
      END IF

      DEALLOCATE (eig)

!      CALL check_diag(matrix, eigenvectors, n)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_dlaf

!***************************************************************************************************
!> \brief DLA-Future eigensolver
!> \author Rocco Meli
!***************************************************************************************************
   SUBROUTINE cp_fm_diag_dlaf_base(matrix, eigenvectors, eigenvalues)
      TYPE(cp_fm_type), INTENT(IN) :: matrix
      TYPE(cp_fm_type), INTENT(IN) :: eigenvectors
      REAL(kind=dp), DIMENSION(:), INTENT(OUT), TARGET :: eigenvalues

      INTEGER :: handle, n, nmo
      INTEGER, TARGET :: info

      INTEGER, DIMENSION(9)                    :: desca, descz
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, z

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_dlaf_base'
      CHARACTER(LEN=100)   :: message

      INTERFACE
         SUBROUTINE pdsyevd_dlaf(n, a, ia, ja, desca, w, z, iz, jz, descz, info) &
            BIND(C, name='dlaf_pdsyevd_wrapper')

            IMPORT :: C_INT, C_DOUBLE, C_CHAR, C_PTR

            INTEGER(kind=C_INT), VALUE :: n, ia, ja, iz, jz
            TYPE(C_PTR), VALUE :: a, w, z
            INTEGER(kind=C_INT), DIMENSION(9) :: desca, descz
            TYPE(C_PTR), VALUE :: info
         END SUBROUTINE pdsyevd_dlaf
      END INTERFACE

      call timeset(routineN, handle)

#if defined(__DLAF)
      n = matrix%matrix_struct%nrow_global

      ! DLAF needs the lower triangular part
      ! Use eigenvectors matrix as workspace
      CALL cp_fm_upper_to_full(matrix, eigenvectors)

      a => matrix%local_data
      z => eigenvectors%local_data

      desca(:) = matrix%matrix_struct%descriptor(:)
      descz(:) = eigenvectors%matrix_struct%descriptor(:)

      info = -1

      CALL pdsyevd_dlaf( &
         n=n, &
         a=C_LOC(a(1, 1)), ia=1, ja=1, desca=desca, &
         w=C_LOC(eigenvalues(1)), &
         z=C_LOC(z(1, 1)), iz=1, jz=1, descz=descz, &
         info=C_LOC(info) &
         )

      IF (info /= 0) THEN
         WRITE (message, "(A,I0,A)") "ERROR in DLAF_PDSYEVD: Eigensolver failed (INFO = ", info, ")"
         CPABORT(TRIM(message))
      END IF
#else
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      MARK_USED(n)
      MARK_USED(nmo)
      MARK_USED(info)
      CPABORT("CP2K compiled without DLAF library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_dlaf_base

END MODULE cp_fm_dlaf_api
