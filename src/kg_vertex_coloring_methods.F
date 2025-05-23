!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits
!>        unsing a vertex coloring algorithm (DSATUR) to find non-interating
!>        subsets, such that two molecules within the same subset have
!>        small/zero overlap (in other words: this molecular pair is not present in
!>        the neighborlist sab_orb for the current value of EPS_DEFAULT)
!> \par History
!>       2012.07 created [Martin Haeufel]
!>         2013.11 Added pair switching and revised Dsatur [Samuel Andermatt]
!> \author Martin Haeufel
! **************************************************************************************************
MODULE kg_vertex_coloring_methods
   USE bibliography,                    ONLY: Andermatt2016,&
                                              Brelaz1979,&
                                              cite_reference
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_min_heap,                     ONLY: cp_heap_fill,&
                                              cp_heap_new,&
                                              cp_heap_pop,&
                                              cp_heap_release,&
                                              cp_heap_reset_node,&
                                              cp_heap_type,&
                                              valt
   USE input_constants,                 ONLY: kg_color_dsatur,&
                                              kg_color_greedy
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kinds,                           ONLY: int_4,&
                                              int_8
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE vertex
      INTEGER :: id = -1
      INTEGER :: color = -1
      INTEGER :: degree = -1 ! degree (number of neighbors)
      INTEGER :: dsat = -1 ! degree of saturation
      LOGICAL, ALLOCATABLE, DIMENSION(:)       :: color_present ! the colors present on neighbors
      TYPE(vertex_p_type), DIMENSION(:), POINTER :: neighbors => NULL()
   END TYPE vertex

   TYPE vertex_p_type
      TYPE(vertex), POINTER :: vertex => NULL()
   END TYPE vertex_p_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_vertex_coloring_methods'

   PUBLIC :: kg_vertex_coloring

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
!> \param pairs ...
!> \param graph ...
! **************************************************************************************************
   SUBROUTINE kg_create_graph(kg_env, pairs, graph)
      TYPE(kg_environment_type), POINTER                 :: kg_env
      INTEGER(KIND=int_4), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(IN)                     :: pairs
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph

      CHARACTER(len=*), PARAMETER                        :: routineN = 'kg_create_graph'

      INTEGER                                            :: handle, i, imol, ineighbor, jmol, kmol, &
                                                            maxdegree, nneighbors, nnodes

      CALL timeset(routineN, handle)

      CALL cite_reference(Andermatt2016)

! The array pairs contains all interacting (overlapping) pairs of molecules.
! It is assumed to be ordered in the following way:
! (1,2), (1,3), (1,4), ..., (2,3), (2,4), ...
! There are no entries (i,i)
! get the number of nodes = total number of molecules

      nnodes = SIZE(kg_env%molecule_set)

      NULLIFY (graph)
      ALLOCATE (graph(nnodes))

      ! allocate and initialize all vertices
      DO i = 1, nnodes
         ALLOCATE (graph(i)%vertex)
         graph(i)%vertex%id = i ! id = imol (molecular index)
         graph(i)%vertex%color = 0 ! means vertex is not colored yet
         graph(i)%vertex%dsat = 0 ! no colored neighbors yet
         graph(i)%vertex%degree = 0 ! as needed for maxdegree....
      END DO

      ! allocate the neighbor lists
      imol = 0

      maxdegree = 0

      DO i = 1, SIZE(pairs, 2)
         jmol = pairs(1, i)
         ! counting loop
         IF (jmol .NE. imol) THEN
            IF (imol .NE. 0) THEN
               ALLOCATE (graph(imol)%vertex%neighbors(nneighbors))
               graph(imol)%vertex%degree = nneighbors
               IF (nneighbors .GT. maxdegree) maxdegree = nneighbors
            END IF
            imol = jmol
            nneighbors = 0
         END IF
         nneighbors = nneighbors + 1
      END DO

      IF (imol .NE. 0) THEN
         ALLOCATE (graph(imol)%vertex%neighbors(nneighbors))
         graph(imol)%vertex%degree = nneighbors
         IF (nneighbors .GT. maxdegree) maxdegree = nneighbors
      END IF

      kg_env%maxdegree = maxdegree

      ! there can be now some nodes that have no neighbors, thus vertex%neighbors
      ! is NOT ASSOCIATED

      ! now add the neighbors - if there are any
      imol = 0
      ineighbor = 0

      DO i = 1, SIZE(pairs, 2)
         jmol = pairs(1, i)
         IF (jmol .NE. imol) THEN
            ineighbor = 0
            imol = jmol
         END IF
         ineighbor = ineighbor + 1
         kmol = pairs(2, i)
         graph(imol)%vertex%neighbors(ineighbor)%vertex => graph(kmol)%vertex
      END DO

      DO i = 1, SIZE(graph)
         IF (graph(i)%vertex%degree > 0) THEN
            ALLOCATE (graph(i)%vertex%color_present(100))
            graph(i)%vertex%color_present(:) = .FALSE.
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE

   ! greedy algorithm
! **************************************************************************************************
!> \brief ...
!> \param graph ...
!> \param maxdegree ...
!> \param ncolors ...
! **************************************************************************************************
   SUBROUTINE color_graph_greedy(graph, maxdegree, ncolors)
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph
      INTEGER, INTENT(in)                                :: maxdegree
      INTEGER, INTENT(out)                               :: ncolors

      CHARACTER(len=*), PARAMETER :: routineN = 'color_graph_greedy'

      INTEGER                                            :: color, handle, i, j, newcolor, &
                                                            nneighbors, nnodes
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: color_present

      CALL timeset(routineN, handle)

      ncolors = 0

      nnodes = SIZE(graph)

      ALLOCATE (color_present(maxdegree + 1))

      DO i = 1, nnodes
         color_present(:) = .FALSE.
         IF (ASSOCIATED(graph(i)%vertex%neighbors)) THEN
            nneighbors = SIZE(graph(i)%vertex%neighbors)
            DO j = 1, nneighbors
               color = graph(i)%vertex%neighbors(j)%vertex%color
               IF (color .NE. 0) color_present(color) = .TRUE.
            END DO
         END IF
         DO j = 1, maxdegree + 1 !nnodes
            IF (color_present(j) .EQV. .FALSE.) THEN
               newcolor = j
               EXIT
            END IF
         END DO
         IF (newcolor .GT. ncolors) ncolors = newcolor
         graph(i)%vertex%color = newcolor ! smallest possible
      END DO

      DEALLOCATE (color_present)

      CALL timestop(handle)

   END SUBROUTINE

   ! prints the subset info to the screen - useful for vmd visualization
   ! note that the index starts with '0' and not with '1'
! **************************************************************************************************
!> \brief ...
!> \param graph ...
!> \param ncolors ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_subsets(graph, ncolors, unit_nr)
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph
      INTEGER, INTENT(IN)                                :: ncolors, unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'print_subsets'

      INTEGER                                            :: counter, handle, i, j, nnodes

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A,T10,A)') "Color |", "Molecules in this subset (IDs start from 0)"

         nnodes = SIZE(graph)

         DO i = 1, ncolors
            WRITE (unit_nr, '(T2,I5,1X,A)', ADVANCE='NO') i, "|"
            counter = 0
            DO j = 1, nnodes
               IF (graph(j)%vertex%color .EQ. i) THEN
                  counter = counter + 1
                  IF (MOD(counter, 13) .EQ. 0) THEN
                     counter = counter + 1
                     WRITE (unit_nr, '()') ! line break
                     WRITE (unit_nr, '(6X,A)', ADVANCE='NO') " |" ! indent next line
                  END IF
                  WRITE (unit_nr, '(I5,1X)', ADVANCE='NO') graph(j)%vertex%id - 1
               END IF
            END DO
            WRITE (unit_nr, '()')
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dsat ...
!> \param degree ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION kg_get_value(dsat, degree) RESULT(value)
      INTEGER, INTENT(IN)                                :: dsat, degree
      INTEGER(KIND=valt)                                 :: value

      INTEGER, PARAMETER                                 :: huge_4 = 2_int_4**30

!   INTEGER, PARAMETER                       :: huge_4 = HUGE(0_int_4) ! PR67219 workaround
! we actually need a max_heap

      value = (huge_4 - INT(dsat, KIND=int_8))*huge_4 + huge_4 - degree

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param graph ...
! **************************************************************************************************
   SUBROUTINE kg_cp_heap_fill(heap, graph)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      TYPE(vertex_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: graph

      INTEGER                                            :: i, nnodes
      INTEGER(kind=valt), ALLOCATABLE, DIMENSION(:)      :: values

      nnodes = SIZE(graph)

      ALLOCATE (values(nnodes))

      DO i = 1, nnodes
         values(i) = kg_get_value(0, graph(i)%vertex%degree)
      END DO

      ! fill the heap
      CALL cp_heap_fill(heap, values)

      DEALLOCATE (values)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param node ...
! **************************************************************************************************
   SUBROUTINE kg_update_node(heap, node)
      TYPE(cp_heap_type)                                 :: heap
      TYPE(vertex), INTENT(IN), POINTER                  :: node

      INTEGER                                            :: degree, dsat, id
      INTEGER(KIND=valt)                                 :: value

      id = node%id

      ! only update node if not yet colored
      IF (heap%index(id) > 0) THEN

         degree = node%degree
         dsat = node%dsat

         value = kg_get_value(dsat, degree)

         CALL cp_heap_reset_node(heap, id, value)

      END IF

   END SUBROUTINE

   ! Subroutine revised by Samuel Andermatt (2013.11)
! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
!> \param graph ...
!> \param ncolors ...
! **************************************************************************************************
   SUBROUTINE kg_dsatur(kg_env, graph, ncolors)
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph
      INTEGER(KIND=int_4), INTENT(OUT)                   :: ncolors

      CHARACTER(len=*), PARAMETER                        :: routineN = 'kg_dsatur'

      INTEGER                                            :: color_limit, handle, i, j, key, &
                                                            maxdegree, nneighbors, nnodes
      INTEGER(KIND=valt)                                 :: value
      LOGICAL                                            :: found
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: color_present
      TYPE(cp_heap_type)                                 :: heap
      TYPE(vertex), POINTER                              :: neighbor, this

      CALL timeset(routineN, handle)

      CALL cite_reference(Brelaz1979)

      ncolors = 0
      color_limit = 100
      maxdegree = kg_env%maxdegree
      nnodes = SIZE(graph)

      IF (kg_env%maxdegree .EQ. 0) THEN
         ! all nodes are disconnected

         ncolors = 1

         DO i = 1, nnodes
            ! only one color needed to color the graph
            graph(i)%vertex%color = 1
         END DO

      ELSE

         ! allocate and fill the heap
         CALL cp_heap_new(heap, nnodes)

         CALL kg_cp_heap_fill(heap, graph)

         DO WHILE (heap%n > 0)

            CALL cp_heap_pop(heap, key, value, found)

            this => graph(key)%vertex

            nneighbors = 0

            IF (ASSOCIATED(this%neighbors)) THEN
               nneighbors = SIZE(this%neighbors)
            ELSE !node is isolated
               this%color = 1
               CYCLE
            END IF
            DO i = 1, this%degree + 1
               IF (this%color_present(i) .EQV. .FALSE.) THEN
                  this%color = i ! smallest possible
                  EXIT
               END IF
            END DO
            IF (this%color .GT. ncolors) ncolors = this%color

            ! update all neighboring nodes
            DO i = 1, nneighbors
               neighbor => this%neighbors(i)%vertex
               IF (neighbor%color_present(this%color)) CYCLE
               neighbor%color_present(this%color) = .TRUE.
               neighbor%dsat = neighbor%dsat + 1
               IF (neighbor%color /= 0) CYCLE
               CALL kg_update_node(heap, neighbor)

            END DO

            IF (this%color == color_limit) THEN !resize all color_present arrays
               ALLOCATE (color_present(color_limit))
               DO i = 1, nnodes
                  IF (graph(i)%vertex%degree == 0) CYCLE
                  color_present(:) = graph(i)%vertex%color_present
                  DEALLOCATE (graph(i)%vertex%color_present)
                  ALLOCATE (graph(i)%vertex%color_present(color_limit*2))
                  DO j = 1, color_limit
                     graph(i)%vertex%color_present(j) = color_present(j)
                  END DO
                  DO j = color_limit + 1, 2*color_limit
                     graph(i)%vertex%color_present(j) = .FALSE.
                  END DO
               END DO
               DEALLOCATE (color_present)
               color_limit = color_limit*2
            END IF

         END DO

         ! release the heap
         CALL cp_heap_release(heap)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Checks if the color of two nodes can be exchanged legally
!> \param this ...
!> \param neighbor ...
!> \param low_col ...
!> \param switchable ...
!> \param ncolors ...
!> \param color_present ...
!> \par History
!>       2013.11 created [Samuel Andermatt]
!> \author Samuel Andermatt
! **************************************************************************************************
   SUBROUTINE kg_check_switch(this, neighbor, low_col, switchable, ncolors, color_present)

      TYPE(vertex), POINTER                              :: this, neighbor
      INTEGER, INTENT(OUT)                               :: low_col
      LOGICAL                                            :: switchable
      INTEGER                                            :: ncolors
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: color_present

      INTEGER                                            :: i

      switchable = .TRUE.
      low_col = ncolors + 1

      DO i = 1, this%degree
         IF (this%neighbors(i)%vertex%id == neighbor%id) CYCLE
         IF (this%neighbors(i)%vertex%color == neighbor%color) THEN
            switchable = .FALSE.
            EXIT
         END IF
      END DO
      IF (.NOT. switchable) RETURN

      color_present(:) = .FALSE.
      color_present(neighbor%color) = .TRUE.
      DO i = 1, neighbor%degree
         IF (neighbor%neighbors(i)%vertex%id == this%id) CYCLE
         color_present(neighbor%neighbors(i)%vertex%color) = .TRUE.
      END DO
      DO i = 1, this%color
         IF (.NOT. color_present(i)) THEN
            low_col = i
            EXIT
         END IF
      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief An algorithm that works on an already colored graph and tries to
!>          reduce the amount of colors by switching the colors of
!>          neighboring vertices
!> \param graph ...
!> \param ncolors ...
!> \par History
!>       2013.11 created [Samuel Andermatt]
!> \author Samuel Andermatt
! **************************************************************************************************
   SUBROUTINE kg_pair_switching(graph, ncolors)
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph
      INTEGER                                            :: ncolors

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_pair_switching'

      INTEGER                                            :: depth, handle, i, iter, j, low_col, &
                                                            low_col_neigh, maxdepth, &
                                                            maxiterations, nnodes, partner
      LOGICAL                                            :: switchable
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: color_present
      TYPE(vertex), POINTER                              :: neighbor, this

      CALL timeset(routineN, handle)
      nnodes = SIZE(graph)
      ALLOCATE (color_present(ncolors))

      !Some default values that should work well
      maxdepth = 4
      maxiterations = 20

      DO iter = 1, maxiterations
         DO depth = 1, maxdepth !First the nodes with larges colors are attempted to be switched and reduced
            !Go trough all nodes and try if you can reduce their color by switching with a neighbour
            DO i = 1, nnodes
               this => graph(i)%vertex
               IF (.NOT. ASSOCIATED(this%neighbors)) CYCLE
               IF (graph(i)%vertex%color < ncolors - depth + 1) CYCLE !Node already has a low color

               partner = 0
               low_col = this%color + 1

               DO j = 1, this%degree
                  neighbor => this%neighbors(j)%vertex
                  IF (neighbor%color > this%color) CYCLE
                  CALL kg_check_switch(this, neighbor, low_col_neigh, switchable, ncolors, color_present)
                  IF (switchable .AND. low_col_neigh < low_col) THEN
                     partner = j
                     low_col = low_col_neigh
                  END IF
               END DO
               IF (partner == 0) CYCLE !Cannot switch favourably with anybody
               !Switch the nodes
               this%color = this%neighbors(partner)%vertex%color
               this%neighbors(partner)%vertex%color = low_col
            END DO

            ncolors = 0
            DO j = 1, nnodes
               IF (graph(j)%vertex%color > ncolors) THEN
                  ncolors = graph(j)%vertex%color
               END IF
            END DO
         END DO
      END DO

      DEALLOCATE (color_present)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param graph ...
!> \param valid ...
! **************************************************************************************************
   SUBROUTINE check_coloring(graph, valid)
      TYPE(vertex_p_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: graph
      LOGICAL, INTENT(out)                               :: valid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'check_coloring'

      INTEGER                                            :: handle, i, j, nneighbors, nnodes
      TYPE(vertex), POINTER                              :: neighbor, node

      CALL timeset(routineN, handle)

      valid = .TRUE.
      nnodes = SIZE(graph)

      DO i = 1, nnodes
         node => graph(i)%vertex
         IF (ASSOCIATED(node%neighbors)) THEN
            nneighbors = SIZE(node%neighbors)
            DO j = 1, nneighbors
               neighbor => node%neighbors(j)%vertex
               IF (neighbor%color == node%color) THEN
                  valid = .FALSE.
                  RETURN
               END IF
            END DO
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param graph ...
! **************************************************************************************************
   SUBROUTINE deallocate_graph(graph)
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph

      INTEGER                                            :: i, nnodes

      nnodes = SIZE(graph)

      DO i = 1, nnodes
         IF (ASSOCIATED(graph(i)%vertex%neighbors)) DEALLOCATE (graph(i)%vertex%neighbors)
         DEALLOCATE (graph(i)%vertex)
      END DO
      DEALLOCATE (graph)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
!> \param pairs ...
!> \param ncolors ...
!> \param color_of_node ...
! **************************************************************************************************
   SUBROUTINE kg_vertex_coloring(kg_env, pairs, ncolors, color_of_node)
      TYPE(kg_environment_type), POINTER                 :: kg_env
      INTEGER(KIND=int_4), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(IN)                     :: pairs
      INTEGER(KIND=int_4), INTENT(OUT)                   :: ncolors
      INTEGER(KIND=int_4), ALLOCATABLE, DIMENSION(:), &
         INTENT(out)                                     :: color_of_node

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_vertex_coloring'

      INTEGER                                            :: color, handle, i, nnodes, unit_nr
      LOGICAL                                            :: valid
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(vertex_p_type), DIMENSION(:), POINTER         :: graph

! get a useful output_unit

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      CALL timeset(routineN, handle)

      CALL kg_create_graph(kg_env, pairs, graph)

      SELECT CASE (kg_env%coloring_method)
      CASE (kg_color_greedy)
         ! simple greedy algorithm
         CALL color_graph_greedy(graph, kg_env%maxdegree, ncolors)
      CASE (kg_color_dsatur)
         ! color with max degree of saturation
         CALL kg_dsatur(kg_env, graph, ncolors)
      CASE DEFAULT
         CPABORT("Coloring method not known.")
      END SELECT

      CALL kg_pair_switching(graph, ncolors)

      valid = .FALSE.
      CALL check_coloring(graph, valid)
      IF (.NOT. valid) &
         CPABORT("Coloring not valid.")

      nnodes = SIZE(kg_env%molecule_set)

      ALLOCATE (color_of_node(nnodes))

      ! gather the subset info in a simple integer array
      DO i = 1, nnodes
         color = graph(i)%vertex%color
         color_of_node(i) = color
      END DO

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 30), " KG coloring result ", REPEAT("-", 29)
         IF (.FALSE.) THEN ! more output for VMD
            WRITE (unit_nr, '()')
            CALL print_subsets(graph, ncolors, unit_nr)
            WRITE (unit_nr, '()')
         END IF
         WRITE (unit_nr, '(T2, A16,50X,I6,1X,A6)') 'Number of colors', ncolors, 'colors'
         IF (valid) WRITE (unit_nr, '(T2, A17,59X,A3)') 'Coloring verified', 'yes'
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

      END IF

      CALL deallocate_graph(graph)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE
