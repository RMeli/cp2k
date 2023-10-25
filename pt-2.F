MODULE cp_dlaf_utils
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#include "src/base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dlaf_utils'

   PUBLIC:: cp_dlaf_initialize

CONTAINS

   SUBROUTINE cp_dlaf_initialzie()
      INTERFACE
         SUBROUTINE dlaf_init_aux() &
            BIND(C, name='dlaf_init')
         END SUBROUTINE dlaf_init_aux
      END INTERFACE

      call dlaf_init_aux()
   END SUBROUTINE cp_dlaf_initialzie
END MODULE cp_dlaf_utils
