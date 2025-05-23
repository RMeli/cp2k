!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief All kind of helpful little routines
!> \par History
!>      none
!> \author CJM & JGH
! **************************************************************************************************
MODULE util
   USE cp_array_sort,                   ONLY: cp_1d_i4_sort,&
                                              cp_1d_i8_sort,&
                                              cp_1d_r_sort,&
                                              cp_1d_s_sort
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'util'
   PUBLIC :: sort, &
             get_limit, &
             locate, &
             find_boundary, &
             sort_unique

   INTERFACE sort
      MODULE PROCEDURE cp_1d_s_sort, cp_1d_r_sort, cp_1d_i4_sort, cp_1d_i8_sort, &
         sort_cv, sort_im, sort_cm
   END INTERFACE

   INTERFACE sort_unique
      MODULE PROCEDURE sort_unique1
   END INTERFACE

   INTERFACE find_boundary
      MODULE PROCEDURE find_boundary1, find_boundary2
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Purpose: Given an array array(1:n), and given a value x, a value x_index
!>             is returned which is the index value of the array element equal
!>             to the value x: x = array(x_index)
!>             The array must be monotonic increasing.
!>             x_index = 0 is returned, if no array element equal to the value
!>             of x was found.
!> \param array ...
!> \param x ...
!> \return ...
!> \par History
!>      Derived from the locate function described in
!>      Numerical Recipes in Fortran 90 (09.01.2004,MK)
! **************************************************************************************************
   PURE FUNCTION locate(array, x) RESULT(x_index)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: array
      INTEGER, INTENT(IN)                                :: x
      INTEGER                                            :: x_index

      INTEGER                                            :: jl, jm, ju, n

      x_index = 0

      IF (x < array(1)) RETURN
      n = SIZE(array)
      IF (x > array(n)) RETURN
      jl = 0
      ju = n + 1
      DO WHILE (ju - jl > 1)
         jm = (ju + jl)/2
         IF (x >= array(jm)) THEN
            jl = jm
         ELSE
            ju = jm
         END IF
      END DO
      IF (x == array(jl)) x_index = jl
   END FUNCTION locate

! **************************************************************************************************
!> \brief Sorts and returns a logical that checks if all elements are unique
!> \param arr ...
!> \param unique ...
!> \par History
!>      Teodoro Laino - Zurich University [tlaino] 04.2007
! **************************************************************************************************
   SUBROUTINE sort_unique1(arr, unique)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: arr
      LOGICAL, INTENT(OUT)                               :: unique

      INTEGER                                            :: i, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: wrk

      n = SIZE(arr)
      unique = .TRUE.
      ALLOCATE (wrk(n))
      CALL sort(arr, n, wrk)
      DO i = 2, n
         IF (arr(i) == arr(i - 1)) THEN
            unique = .FALSE.
            EXIT
         END IF
      END DO
      DEALLOCATE (wrk)
   END SUBROUTINE sort_unique1

! **************************************************************************************************
!> \brief Sorts an array of strings
!> \param arr ...
!> \param n ...
!> \param index ...
!> \author Teodoro Laino [tlaino] - University of Zurich  10.2008
! **************************************************************************************************
   SUBROUTINE sort_cv(arr, n, index)
      INTEGER, INTENT(IN)                                :: n
      CHARACTER(LEN=*), INTENT(INOUT)                    :: arr(1:n)
      INTEGER, INTENT(OUT)                               :: INDEX(1:n)

      INTEGER                                            :: i, j, max_length
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: entries

      max_length = 0
      DO i = 1, n
         max_length = MAX(max_length, LEN_TRIM(arr(i)))
      END DO
      ALLOCATE (entries(max_length, SIZE(arr)))
      DO i = 1, n
         DO j = 1, LEN_TRIM(arr(i))
            entries(j, i) = ICHAR(arr(i) (j:j))
         END DO
         IF (j <= max_length) THEN
            entries(j:max_length, i) = ICHAR(" ")
         END IF
      END DO
      CALL sort_im(entries, istart=1, iend=n, j=1, jsize=max_length, INDEX=INDEX)
      ! Recover string once ordered
      DO i = 1, n
         DO j = 1, max_length
            arr(i) (j:j) = CHAR(entries(j, i))
         END DO
      END DO
      DEALLOCATE (entries)
   END SUBROUTINE sort_cv

#:for argtype, itemtype in [['INTEGER', 'INTEGER'], ['CHARACTER(LEN=*)', 'CHARACTER(LEN=LEN(matrix))']]
! **************************************************************************************************
!> \brief Sorts a multiple arrays M(j,i), ordering iteratively over i with fixed j
!> \param matrix ...
!> \param istart ...
!> \param iend ...
!> \param j ...
!> \param jsize ...
!> \param index ...
!> \author Teodoro Laino [tlaino] - University of Zurich  10.2008
! **************************************************************************************************
   RECURSIVE SUBROUTINE sort_${argtype[0]}$m(matrix, istart, iend, j, jsize, index)
      ${argtype}$, DIMENSION(:, :), INTENT(INOUT)            :: matrix
      INTEGER, INTENT(IN)                                    :: istart, iend, j, jsize
      INTEGER, DIMENSION(:), INTENT(INOUT)                   :: index


      ${itemtype}$                                           :: item
      ${itemtype}$, ALLOCATABLE, DIMENSION(:)                :: work, work2
      INTEGER                                                :: i, ind, isize, k, kend, kstart
      INTEGER, ALLOCATABLE, DIMENSION(:)                     :: bck_index, tmp_index

    isize = iend - istart + 1
    ! Initialize the INDEX array only for the first row..
    IF (j == 1) THEN
       DO i = 1, isize
          INDEX(i) = i
       ENDDO
    END IF

    ! Allocate scratch arrays
    ALLOCATE (work(isize), work2(isize), tmp_index(isize), bck_index(isize))
    ind = 0
    DO i = istart, iend
       ind = ind + 1
       work(ind) = matrix(j, i)
       bck_index(ind) = INDEX(i)
    END DO

    ! Ordering row (j) interval istart..iend
    CALL sort(work, isize, tmp_index)

    ! Copy into global INDEX array with a proper mapping
    ind = 0
    DO i = istart, iend
       ind = ind + 1
       INDEX(i) = bck_index(tmp_index(ind))
       matrix(j, i) = work(ind)
    END DO

    ! Reorder the rest of the array according the present reordering
    DO k = j + 1, jsize
       ind = 0
       DO i = istart, iend
          ind = ind + 1
          work2(ind) = matrix(k, i)
       END DO
       ind = 0
       DO i = istart, iend
          ind = ind + 1
          matrix(k, i) = work2(tmp_index(ind))
       END DO
    END DO

    ! There are more rows to order..
    IF (j < jsize) THEN
       kstart = istart
       item = work(1)
       ind = 0
       DO i = istart, iend
          ind = ind + 1
          IF (item /= work(ind)) THEN
             kend = i - 1
             IF (kstart /= kend) THEN
                CALL sort(matrix, kstart, kend, j + 1, jsize, INDEX)
             END IF
             item = work(ind)
             kstart = i
          END IF
       END DO
       kend = i - 1
       IF (kstart /= kend) THEN
          CALL sort(matrix, kstart, kend, j + 1, jsize, INDEX)
       END IF
    END IF
    DEALLOCATE (work, work2, tmp_index, bck_index)

   END SUBROUTINE sort_${argtype[0]}$m
#:endfor

! **************************************************************************************************
!> \brief divide m entries into n parts, return size of part me
!> \param m ...
!> \param n ...
!> \param me ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_limit(m, n, me) RESULT(nlim)
      INTEGER, INTENT(IN)                                :: m, n, me
      INTEGER                                            :: nlim(2)

      INTEGER                                            :: nl, nu
      REAL(KIND=dp)                                      :: part

      part = REAL(m, KIND=dp)/REAL(n, KIND=dp)
      nl = NINT(REAL(me, KIND=dp)*part) + 1
      nu = NINT(REAL(me + 1, KIND=dp)*part)
      nlim(1) = MAX(1, nl)
      nlim(2) = MIN(m, nu)

   END FUNCTION get_limit

! **************************************************************************************************
!> \brief finds boundary where element search starts and ends in a 1D array
!>      array1:      XXXXXAAAAAAAAAXXDGFSFGWDDDDDDDAAAWE
!>                        |       |
!>                     start     end  (searching for A)
!> \param num_array ...
!> \param ntot ...
!> \param first ...
!> \param last ...
!> \param search ...
! **************************************************************************************************
   PURE SUBROUTINE find_boundary1(num_array, ntot, first, last, search)
      INTEGER, POINTER                                   :: num_array(:)
      INTEGER, INTENT(IN)                                :: ntot
      INTEGER, INTENT(OUT)                               :: first, last
      INTEGER, INTENT(IN)                                :: search

      INTEGER                                            :: i
      LOGICAL                                            :: found

      found = .FALSE.
      first = 0
      last = ntot

      DO i = 1, ntot
         IF (num_array(i) == search) THEN
            IF (.NOT. found) THEN
               first = i
            END IF
            found = .TRUE.
         ELSE
            IF (found) THEN
               last = i - 1
               EXIT
            END IF
            found = .FALSE.
         END IF
      END DO

   END SUBROUTINE find_boundary1

! **************************************************************************************************
!> \brief finds boundary where element search1 starts and ends in array1 checking
!>      at the same time search2 in array2
!>      array1:      XXXXXAAAAAAAAAXXDGFSFGWDDDDDDDAAAWE
!>      array2:      XXXXASDEYYYYASDEFAAAARGASGASRGAWRRR
!>                           |  |
!>                       start  end  (searching for A and Y)
!> \param num_array1 ...
!> \param num_array2 ...
!> \param ntot ...
!> \param first ...
!> \param last ...
!> \param search1 ...
!> \param search2 ...
! **************************************************************************************************
   PURE SUBROUTINE find_boundary2(num_array1, num_array2, ntot, first, last, search1, search2)
      INTEGER, POINTER                                   :: num_array1(:), num_array2(:)
      INTEGER, INTENT(IN)                                :: ntot
      INTEGER, INTENT(OUT)                               :: first, last
      INTEGER, INTENT(IN)                                :: search1, search2

      INTEGER                                            :: i, tfirst, tlast
      LOGICAL                                            :: found

      found = .FALSE.
      first = 0
      last = ntot

      CALL find_boundary(num_array1, ntot, tfirst, tlast, search1)
      last = tlast
      DO i = tfirst, tlast
         IF (num_array2(i) == search2) THEN
            IF (.NOT. found) THEN
               first = i
            END IF
            found = .TRUE.
         ELSE
            IF (found) THEN
               last = i - 1
               EXIT
            END IF
            found = .FALSE.
         END IF
      END DO

   END SUBROUTINE find_boundary2

END MODULE util
