MODULE foo
   USE ISO_C_BINDING, ONLY: C_INT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: bar

CONTAINS

   SUBROUTINE bar()
      INTEGER :: i

      INTERFACE
         SUBROUTINE bar_inter(ii) BIND(C, name="bar_c")
            IMPORT C_INT
            INTEGER(kind=C_INT), VALUE :: ii
         END SUBROUTINE bar_inter
      END INTERFACE

      CALL bar_inter(i)

   END SUBROUTINE bar
END MODULE foo
