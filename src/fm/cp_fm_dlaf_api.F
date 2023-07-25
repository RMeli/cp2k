MODULE cp_fm_dlaf_api

   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_CHAR, C_FLOAT, C_LOC

   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_types, ONLY: cp_fm_type
   USE kinds, ONLY: dp

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_dlaf_api'

   PUBLIC :: cp_dlaf_init, cp_dlaf_finalize
   PUBLIC :: cp_pdpotrf_dlaf, cp_pspotrf_dlaf
   !PUBLIC :: cp_pdsyevd_dlaf
   PUBLIC :: cp_fm_diag_dlaf

CONTAINS

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
      INTEGER, dimension(:), TARGET       ::        DESCA
      REAL(KIND=dp), DIMENSION(:, :), TARGET   :: A

      INTERFACE
         SUBROUTINE pdpotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pdpotrf')
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
            CHARACTER(KIND=C_CHAR), VALUE :: UPLO
            INTEGER(KIND=C_INT), VALUE      ::        IA, JA, N
            TYPE(C_PTR), value ::       INFO
            TYPE(C_PTR), value ::        DESCA
            TYPE(C_PTR), value ::              A
         END SUBROUTINE pdpotrf_dlaf
      END INTERFACE

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(A))
#endif
      CALL cp_dlaf_create_grid(desca(2))
      CALL pdpotrf_dlaf(UPLO, N, C_LOC(A(1, 1)), IA, JA, C_LOC(DESCA(1)), C_LOC(INFO))
      CALL cp_dlaf_free_grid(desca(2))
   END SUBROUTINE cp_pdpotrf_dlaf

   SUBROUTINE cp_pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO)
      CHARACTER, INTENT(IN) ::           UPLO
      INTEGER, INTENT(IN)       ::       IA, JA, N
      INTEGER, TARGET     ::      INFO
      INTEGER, dimension(:), TARGET      ::        DESCA
      REAL, dimension(:, :), TARGET    ::              A
      INTERFACE
         SUBROUTINE pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
            BIND(C, name='dlaf_pspotrf')
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
            CHARACTER(KIND=C_CHAR), VALUE :: UPLO
            INTEGER(KIND=C_INT), value      ::        IA, JA, N
            TYPE(C_PTR), value    ::       INFO
            TYPE(C_PTR), value      ::        DESCA
            TYPE(C_PTR), value   ::              A
         END SUBROUTINE pspotrf_dlaf
      END INTERFACE
      CALL cp_dlaf_create_grid(desca(2))
      CALL pspotrf_dlaf(UPLO, N, C_LOC(A(1, 1)), IA, JA, C_LOC(DESCA(1)), C_LOC(INFO))
      CALL cp_dlaf_free_grid(desca(2))
   END SUBROUTINE cp_pspotrf_dlaf

!    ! Compute all eigenvalues and eigenvectors of a real symmetric matrix
!    SUBROUTINE cp_pdsyevd_dlaf(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info)
!       CHARACTER, INTENT(IN) :: jobz
!       CHARACTER, INTENT(IN) :: uplo
!       INTEGER, INTENT(IN) :: n, ia, ja, iz, jz
!       REAL(kind=dp), DIMENSION(:), TARGET :: work
!       INTEGER, DIMENSION(:), TARGET :: iwork
!       INTEGER, INTENT(IN) :: lwork, liwork
!       INTEGER, TARGET :: info
!       INTEGER, DIMENSION(:), TARGET :: desca, descz
!       REAL(KIND=dp), DIMENSION(:), TARGET :: w
!       REAL(KIND=dp), DIMENSION(:, :), TARGET :: a, z
!       INTERFACE
!          SUBROUTINE pdsyevd_dlaf(UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, INFO) &
!             BIND(C, name='dlaf_pdsyevd')
!             IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
!             CHARACTER(KIND=C_CHAR), VALUE :: UPLO
!             INTEGER(KIND=C_INT), VALUE :: N, IA, JA, IZ, JZ
!             TYPE(C_PTR), value :: INFO
!             TYPE(C_PTR), value :: DESCA, DESCZ
!             TYPE(C_PTR), value :: A, W, Z
!          END SUBROUTINE pdsyevd_dlaf
!       END INTERFACE
! #if __GNUC__ >= 9
!       CPASSERT(IS_CONTIGUOUS(A))
! #endif
!       CALL cp_dlaf_create_grid(desca(2))
!       CALL pdsyevd_dlaf( &
!          "L", n, C_LOC(a(1, 1)), ia, ja, C_LOC(desca(1)), &
!          C_LOC(w(1)), C_LOC(z(1, 1)), iz, jz, C_LOC(descz(1)), &
!          C_LOC(info) &
!          )
!       CALL cp_dlaf_free_grid(desca(2))
!    END SUBROUTINE cp_pdsyevd_dlaf

SUBROUTINE cp_fm_diag_dlaf(matrix, eigenvectors, eigenvalues, info)
   TYPE(cp_fm_type), INTENT(IN) :: matrix
   TYPE(cp_fm_type), INTENT(IN) :: eigenvectors
   REAL(kind=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues
   INTEGER, INTENT(OUT) :: info

   REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: eigenvalues_buffer

   INTEGER :: handle, n, nmo, blacs_context

   CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_dlaf'


   INTERFACE
      SUBROUTINE pdsyevd_dlaf(uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, info) &
         BIND(C, name='dlaf_pdsyevd')

         IMPORT :: C_INT, C_DOUBLE, C_CHAR
         CHARACTER(kind=C_CHAR), VALUE :: uplo
         INTEGER(kind=C_INT), VALUE :: n, ia, ja, iz, jz, info
         REAL(kind=C_DOUBLE), DIMENSION(*) :: a, w, z
         INTEGER(kind=C_INT), DIMENSION(*) :: desca, descz
      END SUBROUTINE pdsyevd_dlaf
   END INTERFACE

   call timeset(routineN, handle)
  
#if defined(__DLAF)
   n = matrix%matrix_struct%nrow_global
   blacs_context = matrix%matrix_struct%descriptor(2)

   ! The passed eigenvalues array might be smaller than n
   ALLOCATE(eigenvalues_buffer(n))

   CALL cp_dlaf_create_grid(blacs_context)

   ! TODO: Upper
   CALL pdsyevd_dlaf(&
      uplo="L", &
      n=n, &
      a=matrix%local_data, ia=1, ja=1, desca=matrix%matrix_struct%descriptor, &
      w=eigenvalues_buffer, &
      z=eigenvectors%local_data, iz=1, jz=1, descz=eigenvectors%matrix_struct%descriptor, &
      info=info &
   )
      
   CALL cp_dlaf_free_grid(blacs_context)

   nmo = SIZE(eigenvalues)
   eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)

#else
   MARK_USED(matrix)
   MARK_USED(eigenvectors)
   MARK_USED(eigenvvalues)
   MARK_USED(eigenvvalues_buffer)
   MARK_USED(n)
   MARK_USED(nmo)
   MARK_USED(info)
   MARK_USED(context)
   CPABORT("CP2K compiled without DLAF library.")
#endif

   CALL timestop(handle)
   
END SUBROUTINE cp_fm_diag_dlaf

END MODULE cp_fm_dlaf_api
