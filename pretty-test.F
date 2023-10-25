!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dlaf_utils
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#include "src/base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "cp_dlaf_utils"

   PUBLIC :: cp_dlaf_initialize

CONTAINS

! **************************************************************************************************
!> \brief Initialize DLA-Future and pika runtime
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_initialize()
      INTERFACE
         SUBROUTINE dlaf_init_aux() &
            BIND(C, name="dlaf_init")
         END SUBROUTINE dlaf_init_aux
      END INTERFACE

#if defined(__DLAF)
      call dlaf_init_aux()
#else
      CPABORT("CP2K compiled without the DLA-Future library.")
#endif
   END SUBROUTINE cp_dlaf_initialize

END MODULE cp_dlaf_utils
