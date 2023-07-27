!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_fm_dlaf_api

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_CHAR, C_FLOAT, C_LOC

   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_types, ONLY: cp_fm_type
   USE kinds, ONLY: dp

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_dlaf_api'

#if defined(__DLAF)
   PUBLIC :: cp_dlaf_init, cp_dlaf_finalize
   PUBLIC :: cp_pdpotrf_dlaf, cp_pspotrf_dlaf
#endif
   PUBLIC :: cp_fm_diag_dlaf

CONTAINS
      
#if defined(__DLAF)
   SUBROUTINE cp_dlaf_init()
      INTERFACE
         SUBROUTINE dlaf_init_aux() &
            BIND(C, name='dlaf_init')
         end SUBROUTINE dlaf_init_aux
      end INTERFACE

      call dlaf_init_aux()
   end SUBROUTINE cp_dlaf_init

   SUBROUTINE cp_dlaf_finalize()
      INTERFACE
         SUBROUTINE dlaf_finalize_aux() &
            BIND(C, name='dlaf_finalize')
         end SUBROUTINE dlaf_finalize_aux
      END INTERFACE
      CALL dlaf_finalize_aux()
   END SUBROUTINE cp_dlaf_finalize

   SUBROUTINE cp_dlaf_create_grid(blacs_context)
      INTEGER, INTENT(IN) :: blacs_context
      INTERFACE
         SUBROUTINE dlaf_create_grid(blacs_contxt) &
            BIND(C, name='dlaf_create_grid_from_blacs')
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE :: blacs_contxt
         END SUBROUTINE
      END INTERFACE
      CALL dlaf_create_grid(blacs_context)
   END SUBROUTINE cp_dlaf_create_grid

   SUBROUTINE cp_dlaf_free_grid(blacs_context)
      INTEGER, INTENT(IN) :: blacs_context
      INTERFACE
         SUBROUTINE dlaf_free_grid(blacs_contxt) &
            BIND(C, name='dlaf_free_grid')
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE :: blacs_contxt
         END SUBROUTINE
      END INTERFACE
      CALL dlaf_free_grid(blacs_context)
   end SUBROUTINE cp_dlaf_free_grid

   SUBROUTINE cp_pdpotrf_dlaf(uplo, n, a, ia, ja, desca, info)
      CHARACTER, INTENT(IN) ::           uplo
      INTEGER, INTENT(IN)       ::       ia, ja, n
      INTEGER, TARGET     ::      info
      INTEGER, dimension(9)       ::        DESCA
      REAL(KIND=dp), DIMENSION(:, :), TARGET   :: A
      CHARACTER :: myuplo

      INTERFACE
         SUBROUTINE pdpotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pdpotrf')
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
            CHARACTER(KIND=C_CHAR), VALUE :: UPLO
            INTEGER(KIND=C_INT), VALUE      ::        IA, JA, N
            TYPE(C_PTR), value ::       INFO
            INTEGER(kind=C_INT), DIMENSION(*)         :: DESCA
            TYPE(C_PTR), value ::              A
         END SUBROUTINE pdpotrf_dlaf
      END INTERFACE

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(A))
#endif

      ! Passing CHARACTER, INTENT(IN) directly does not work, the
      ! character is shifted by one
      ! TODO: Figure out what's going on
      myuplo = uplo

      CALL cp_dlaf_create_grid(desca(2))
      CALL pdpotrf_dlaf(myuplo, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      CALL cp_dlaf_free_grid(desca(2))
   END SUBROUTINE cp_pdpotrf_dlaf

   SUBROUTINE cp_pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO)
      CHARACTER, INTENT(IN) ::           UPLO
      INTEGER, INTENT(IN)       ::       IA, JA, N
      INTEGER, TARGET     ::      INFO
      INTEGER, dimension(9)      ::        DESCA
      REAL, dimension(:, :), TARGET    ::              A
      CHARACTER :: myuplo
      INTERFACE
         SUBROUTINE pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pspotrf')
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
            CHARACTER(KIND=C_CHAR), VALUE :: UPLO
            INTEGER(KIND=C_INT), value      ::        IA, JA, N
            TYPE(C_PTR), value    ::       INFO
            INTEGER(kind=C_INT), DIMENSION(*)         :: DESCA
            TYPE(C_PTR), value   ::              A
         END SUBROUTINE pspotrf_dlaf
      END INTERFACE

      ! Passing CHARACTER, INTENT(IN) directly does not work, the
      ! character is shifted by one
      ! TODO: Figure out what's going on
      myuplo = uplo

      CALL cp_dlaf_create_grid(desca(2))
      CALL pspotrf_dlaf(myuplo, N, C_LOC(A(1, 1)), IA, JA, DESCA, C_LOC(INFO))
      CALL cp_dlaf_free_grid(desca(2))
   END SUBROUTINE cp_pspotrf_dlaf

#endif

   SUBROUTINE cp_fm_diag_dlaf(matrix, eigenvectors, eigenvalues, info)
      TYPE(cp_fm_type), INTENT(IN) :: matrix
      TYPE(cp_fm_type), INTENT(IN) :: eigenvectors
      REAL(kind=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
      INTEGER, TARGET :: info

      REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: eigenvalues_buffer

      INTEGER :: handle, n, nmo, blacs_context

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_dlaf'

      INTERFACE
         SUBROUTINE pdsyevd_dlaf(uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, info) &
            BIND(C, name='dlaf_pdsyevd')

            IMPORT :: C_INT, C_DOUBLE, C_CHAR, C_PTR
            CHARACTER(kind=C_CHAR), VALUE :: uplo
            INTEGER(kind=C_INT), VALUE :: n, ia, ja, iz, jz
            REAL(kind=C_DOUBLE), DIMENSION(*) :: a, w, z
            INTEGER(kind=C_INT), DIMENSION(*) :: desca, descz
            TYPE(C_PTR) :: info
         END SUBROUTINE pdsyevd_dlaf
      END INTERFACE

      call timeset(routineN, handle)

#if defined(__DLAF)
      n = matrix%matrix_struct%nrow_global
      blacs_context = matrix%matrix_struct%descriptor(2)

      ! The passed eigenvalues array might be smaller than n
      ALLOCATE (eigenvalues_buffer(n))

      CALL cp_dlaf_create_grid(blacs_context)

      ! TODO: Upper
      CALL pdsyevd_dlaf( &
         uplo="L", &
         n=n, &
         a=matrix%local_data, ia=1, ja=1, desca=matrix%matrix_struct%descriptor, &
         w=eigenvalues_buffer, &
         z=eigenvectors%local_data, iz=1, jz=1, descz=eigenvectors%matrix_struct%descriptor, &
         info=C_LOC(info) &
         )

      CALL cp_dlaf_free_grid(blacs_context)

      nmo = SIZE(eigenvalues)
      eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)

#else
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      MARK_USED(eigenvalues_buffer)
      MARK_USED(n)
      MARK_USED(nmo)
      MARK_USED(info)
      CPABORT("CP2K compiled without DLAF library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_dlaf

END MODULE cp_fm_dlaf_api
