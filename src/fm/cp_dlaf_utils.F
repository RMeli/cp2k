!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dlaf_utils

#include "../base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dlaf_utils'

   PUBLIC :: cp_dlaf_initialize, cp_dlaf_finalize
   PUBLIC :: cp_dlaf_create_grid, cp_dlaf_free_grid

CONTAINS

!***************************************************************************************************
!> \brief Initialize DLA-Future and pika runtime
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_initialize()
      END INTERFACE

#if defined(__DLAF)
      CALL dlaf_init_aux()
#endif
      END SUBROUTINE cp_dlaf_initialize

!***************************************************************************************************
!> \brief Finalize DLA-Future and pika runtime
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
!***************************************************************************************************
      SUBROUTINE cp_dlaf_finalize()
         END INTERFACE

#if defined(__DLAF)
         CALL dlaf_finalize_aux()
#endif
      END SUBROUTINE cp_dlaf_finalize

!***************************************************************************************************
!> \brief Create DLA-Future grid from BLACS context
!> \param blacs_context ...
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
      SUBROUTINE cp_dlaf_create_grid(blacs_context)
      INTEGER, INTENT(IN)                                :: blacs_context

         END INTERFACE

#if defined(__DLAF)
         CALL dlaf_create_grid(blacs_context)
#else
         MARK_USED(blacs_context)
#endif
      END SUBROUTINE cp_dlaf_create_grid

!***************************************************************************************************
!> \brief Free DLA-Future grid corresponding to BLACS context
!> \param blacs_context ...
!> \author Rocco Meli
!> \author Mikeal Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
      SUBROUTINE cp_dlaf_free_grid(blacs_context)
      INTEGER, INTENT(IN)                                :: blacs_context

         END INTERFACE

#if defined(__DLAF)
         CALL dlaf_free_grid(blacs_context)
#else
         MARK_USED(blacs_context)
#endif
      END SUBROUTINE cp_dlaf_free_grid

   END MODULE cp_dlaf_utils
