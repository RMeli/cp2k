!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE manybody_quip

   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE atomic_kind_types, ONLY: atomic_kind_type
   USE bibliography, ONLY: QUIP_ref, &
                           cite_reference
   USE cell_types, ONLY: cell_type
   USE message_passing, ONLY: mp_para_env_type
   USE fist_nonbond_env_types, ONLY: fist_nonbond_env_get, &
                                     fist_nonbond_env_set, &
                                     fist_nonbond_env_type, &
                                     quip_data_type
   USE kinds, ONLY: dp
   USE pair_potential_types, ONLY: pair_potential_pp_type, &
                                   pair_potential_single_type, &
                                   quip_pot_type, &
                                   quip_type
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: angstrom, &
                      evolt
#ifdef __QUIP
   USE quip_unified_wrapper_module, ONLY: quip_unified_wrapper
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC quip_energy_store_force_virial, quip_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_quip'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param fist_nonbond_env ...
!> \param pot_quip ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE quip_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, fist_nonbond_env, &
                                             pot_quip, para_env)
      TYPE(particle_type), POINTER             :: particle_set(:)
      TYPE(cell_type), POINTER                 :: cell
      TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER    :: potparm
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(kind=dp)                            :: pot_quip
      TYPE(mp_para_env_type), OPTIONAL, &
         POINTER                                :: para_env

#ifdef __QUIP
      CHARACTER(len=2), ALLOCATABLE            :: elem_symbol(:)
      INTEGER                                  :: i, iat, iat_use, ikind, &
                                                  jkind, n_atoms, n_atoms_use, &
                                                  output_unit
      LOGICAL                                  :: do_parallel
      LOGICAL, ALLOCATABLE                     :: use_atom(:)
      REAL(kind=dp)                            :: lattice(3, 3), virial(3, 3)
      REAL(kind=dp), ALLOCATABLE               :: force(:, :), pos(:, :)
      TYPE(pair_potential_single_type), &
         POINTER                                :: pot
      TYPE(quip_data_type), POINTER            :: quip_data
      TYPE(quip_pot_type), POINTER             :: quip

#endif
#ifndef __QUIP
      MARK_USED(particle_set)
      MARK_USED(cell)
      MARK_USED(atomic_kind_set)
      MARK_USED(potparm)
      MARK_USED(fist_nonbond_env)
      MARK_USED(pot_quip)
      MARK_USED(para_env)
      CALL cp_abort(__LOCATION__, "In order to use QUIP you need to download "// &
                    "and install the libAtoms/QUIP library (check CP2K manual for details)")
#else
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      use_atom = .FALSE.

      NULLIFY (quip)

      DO ikind = 1, SIZE(atomic_kind_set)
      DO jkind = 1, SIZE(atomic_kind_set)
         pot => potparm%pot(ikind, jkind)%pot
         DO i = 1, SIZE(pot%type)
            IF (pot%type(i) /= quip_type) CYCLE
            IF (.NOT. ASSOCIATED(quip)) quip => pot%set(i)%quip
            DO iat = 1, n_atoms
               IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
                   particle_set(iat)%atomic_kind%kind_number == jkind) use_atom(iat) = .TRUE.
            END DO ! iat
         END DO ! i
      END DO ! jkind
      END DO ! ikind
      n_atoms_use = COUNT(use_atom)
      ALLOCATE (pos(3, n_atoms_use), force(3, n_atoms_use), elem_symbol(n_atoms_use))

      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         pos(1:3, iat_use) = particle_set(iat)%r*angstrom
         elem_symbol(iat_use) = particle_set(iat)%atomic_kind%element_symbol
      END DO
      IF (iat_use > 0) CALL cite_reference(QUIP_ref)
      output_unit = cp_logger_get_default_io_unit()
      lattice = cell%hmat*angstrom
      do_parallel = .FALSE.
      IF (PRESENT(para_env)) THEN
         do_parallel = para_env%num_pe > 1
      END IF
      IF (do_parallel) THEN
         CALL quip_unified_wrapper( &
            N=n_atoms_use, pos=pos, lattice=lattice, symbol=elem_symbol, &
            quip_param_file=TRIM(quip%quip_file_name), &
            quip_param_file_len=LEN_TRIM(quip%quip_file_name), &
            init_args_str=TRIM(quip%init_args), &
            init_args_str_len=LEN_TRIM(quip%init_args), &
            calc_args_str=TRIM(quip%calc_args), &
            calc_args_str_len=LEN_TRIM(quip%calc_args), &
            energy=pot_quip, force=force, virial=virial, &
            output_unit=output_unit, mpi_communicator=para_env%get_handle())
      ELSE
         CALL quip_unified_wrapper( &
            N=n_atoms_use, pos=pos, lattice=lattice, symbol=elem_symbol, &
            quip_param_file=TRIM(quip%quip_file_name), &
            quip_param_file_len=LEN_TRIM(quip%quip_file_name), &
            init_args_str=TRIM(quip%init_args), &
            init_args_str_len=LEN_TRIM(quip%init_args), &
            calc_args_str=TRIM(quip%calc_args), &
            calc_args_str_len=LEN_TRIM(quip%calc_args), &
            energy=pot_quip, force=force, virial=virial, output_unit=output_unit)
      END IF
      ! convert units
      pot_quip = pot_quip/evolt
      force = force/(evolt/angstrom)
      virial = virial/evolt
      ! account for double counting from multiple MPI processes
      IF (PRESENT(para_env)) pot_quip = pot_quip/REAL(para_env%num_pe, dp)
      IF (PRESENT(para_env)) force = force/REAL(para_env%num_pe, dp)
      IF (PRESENT(para_env)) virial = virial/REAL(para_env%num_pe, dp)
      ! get quip_data to save force, virial info
      CALL fist_nonbond_env_get(fist_nonbond_env, quip_data=quip_data)
      IF (.NOT. ASSOCIATED(quip_data)) THEN
         ALLOCATE (quip_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, quip_data=quip_data)
         NULLIFY (quip_data%use_indices, quip_data%force)
      END IF
      IF (ASSOCIATED(quip_data%force)) THEN
         IF (SIZE(quip_data%force, 2) /= n_atoms_use) THEN
            DEALLOCATE (quip_data%force, quip_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(quip_data%force)) THEN
         ALLOCATE (quip_data%force(3, n_atoms_use))
         ALLOCATE (quip_data%use_indices(n_atoms_use))
      END IF
      ! save force, virial info
      iat_use = 0
      DO iat = 1, n_atoms
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            quip_data%use_indices(iat_use) = iat
         END IF
      END DO
      quip_data%force = force
      quip_data%virial = virial

      DEALLOCATE (use_atom, pos, force, elem_symbol)
#endif
   END SUBROUTINE quip_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE quip_add_force_virial(fist_nonbond_env, force, virial)
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(KIND=dp)                            :: force(:, :), virial(3, 3)

#ifdef __QUIP
      INTEGER                                  :: iat, iat_use
      TYPE(quip_data_type), POINTER            :: quip_data
#endif

#ifndef __QUIP
      MARK_USED(fist_nonbond_env)
      MARK_USED(force)
      MARK_USED(virial)
      RETURN
#else
      CALL fist_nonbond_env_get(fist_nonbond_env, quip_data=quip_data)
      IF (.NOT. ASSOCIATED(quip_data)) RETURN

      DO iat_use = 1, SIZE(quip_data%use_indices)
         iat = quip_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(force, 2))
         force(1:3, iat) = force(1:3, iat) + quip_data%force(1:3, iat_use)
      END DO
      virial = virial + quip_data%virial
#endif
   END SUBROUTINE quip_add_force_virial

END MODULE manybody_quip
