MODULE dlaf
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dlaf'
   PUBLIC :: cp_dlaf_init, cp_dlaf_finalize
   PUBLIC :: cp_pdpotrf_dlaf, cp_pspotrf_dlaf
   PUBLIC :: cp_pdsyevd_dlaf
CONTAINS

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE cp_dlaf_init()
      INTERFACE
         SUBROUTINE dlaf_init_aux() &
            BIND(C, name='dlaf_init')
         END SUBROUTINE dlaf_init_aux
      END INTERFACE

      !character(:), allocatable :: arg
      !integer arglen, stat
      !call get_command_argument(number=1, length=arglen)  ! Assume for simplicity success
      !allocate (character(arglen) :: arg)
      !call get_command_argument(number=1, value=arg, status=stat)

      CALL dlaf_init_aux()
   END SUBROUTINE cp_dlaf_init

   SUBROUTINE cp_dlaf_finalize()
      INTERFACE
         SUBROUTINE dlaf_finalize_aux() &
            BIND(C, name='dlaf_finalize')
         END SUBROUTINE dlaf_finalize_aux
      END INTERFACE

      CALL dlaf_finalize_aux()
   END SUBROUTINE cp_dlaf_finalize

   SUBROUTINE cp_pdpotrf_dlaf(uplo, n, a, ia, ja, desca, info)
      CHARACTER, INTENT(IN)                              :: uplo
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), DIMENSION(:, :), TARGET             :: A
      INTEGER, INTENT(IN)                                :: ia, ja
      INTEGER, DIMENSION(:), TARGET                      :: DESCA
      INTEGER, TARGET                                    :: info

      END INTERFACE
#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(A))
#endif
      CALL pdpotrf_dlaf(UPLO, N, C_LOC(A(1, 1)), IA, JA, C_LOC(DESCA(1)), C_LOC(INFO))
      END SUBROUTINE cp_pdpotrf_dlaf

! **************************************************************************************************
!> \brief ...
!> \param UPLO ...
!> \param N ...
!> \param A ...
!> \param IA ...
!> \param JA ...
!> \param DESCA ...
!> \param INFO ...
! **************************************************************************************************
      SUBROUTINE cp_pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO)
      CHARACTER, INTENT(IN)                              :: UPLO
      INTEGER, INTENT(IN)                                :: N
      REAL, DIMENSION(:, :), TARGET                      :: A
      INTEGER, INTENT(IN)                                :: IA, JA
      INTEGER, DIMENSION(:), TARGET                      :: DESCA
      INTEGER, TARGET                                    :: INFO

         INTERFACE
            SUBROUTINE pspotrf_dlaf(UPLO, N, A, IA, JA, DESCA, INFO) &
               BIND(C, name='pspotrf_dlaf')
               IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
               CHARACTER(KIND=C_CHAR), VALUE :: UPLO
               INTEGER(KIND=C_INT), value      ::        IA, JA, N
               TYPE(C_PTR), value    ::       INFO
               TYPE(C_PTR), value      ::        DESCA
               TYPE(C_PTR), value   ::              A
            END SUBROUTINE pspotrf_dlaf
         END INTERFACE

         CALL pspotrf_dlaf(UPLO, N, C_LOC(A(1, 1)), IA, JA, C_LOC(DESCA(1)), C_LOC(INFO))
      END SUBROUTINE cp_pspotrf_dlaf

      ! Compute all eigenvalues and eigenvectors of a real symmetric matrix
! **************************************************************************************************
!> \brief ...
!> \param jobz ...
!> \param uplo ...
!> \param n ...
!> \param a ...
!> \param ia ...
!> \param ja ...
!> \param desca ...
!> \param w ...
!> \param z ...
!> \param iz ...
!> \param jz ...
!> \param descz ...
!> \param work ...
!> \param lwork ...
!> \param iwork ...
!> \param liwork ...
!> \param info ...
! **************************************************************************************************
      SUBROUTINE cp_pdsyevd_dlaf(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info)
      CHARACTER, INTENT(IN)                              :: jobz, uplo
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), DIMENSION(:, :), TARGET             :: a
      INTEGER, INTENT(IN)                                :: ia, ja
      INTEGER, DIMENSION(:), TARGET                      :: desca
      REAL(kind=dp), DIMENSION(:), TARGET                :: w
      REAL(KIND=dp), DIMENSION(:, :), TARGET             :: z
      INTEGER, INTENT(IN)                                :: iz, jz
      INTEGER, DIMENSION(:), TARGET                      :: descz
      REAL(kind=dp), DIMENSION(:), TARGET                :: work
      INTEGER, INTENT(IN)                                :: lwork
      INTEGER, DIMENSION(:), TARGET                      :: iwork
      INTEGER, INTENT(IN)                                :: liwork
      INTEGER, TARGET                                    :: info

      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: printwork

! TODO: Remove
         END INTERFACE
#if __GNUC__ >= 9
         CPASSERT(IS_CONTIGUOUS(A))
#endif
         ! The input matrix is 'U', needst to be transposed to 'L' for DLAF eigensolver
         ALLOCATE (printwork(desca(5)*2))
         ! call pdlaprnt(1, 1, a(1,1), 1, 30, desca(1), 0, 0, "A_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, a(1,1), 30, 1, desca(1), 0, 0, "A_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, z(1,1), 1, 30, descz(1), 0, 0, "Z_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, z(1,1), 30, 1, descz(1), 0, 0, "Z_EIGENSOLVER", 6, printwork)

         ! allocate(aa(n, n), zz(n, n))
         ! call pdtran(n, n, 1.0_dp, a, ia, ja, desca(1), 0.0_dp, aa, ia, ja, desca(1))
         ! call pdtran(n, n, 1.0_dp, z, iz, jz, descz(1), 0.0_dp, zz, iz, jz, descz(1))

         ! call pdlaprnt(1, 1, aa(1,1), 1, 30, desca(1), 0, 0, "AA_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, aa(1,1), 30, 1, desca(1), 0, 0, "AA_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, zz(1,1), 1, 30, descz(1), 0, 0, "ZZ_EIGENSOLVER", 6, printwork)
         ! call pdlaprnt(1, 1, zz(1,1), 30, 1, descz(1), 0, 0, "ZZ_EIGENSOLVER", 6, printwork)

         CALL pdsyevd_dlaf( &
            jobz, uplo, n, C_LOC(a(1, 1)), ia, ja, C_LOC(desca(1)), &
            C_LOC(w(1)), C_LOC(z(1, 1)), iz, jz, C_LOC(descz(1)), &
            C_LOC(work(1)), lwork, C_LOC(iwork(1)), liwork, &
            C_LOC(info) &
            )

         ! CALL pdsyevd_dlaf(&
         !     jobz, uplo, n, C_LOC(aa(1,1)), ia, ja, C_LOC(desca(1)), &
         !     C_LOC(w(1)), C_LOC(zz(1,1)), iz, jz, C_LOC(descz(1)), &
         !     C_LOC(work(1)), lwork, C_LOC(iwork(1)), liwork, &
         !     C_LOC(info) &
         ! )

         ! The output matrix is 'L', needs to be transposed to 'U' for CP2K
         ! call pdtran(n, n, 1.0_dp, zz, iz, jz, descz(1), 0.0_dp, z, iz, jz, descz(1))
         !
         ! deallocate(aa, zz)
      END SUBROUTINE cp_pdsyevd_dlaf

   END MODULE dlaf
