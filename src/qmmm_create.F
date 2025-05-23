!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize a QM/MM calculation
!> \par History
!>      5.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qmmm_create
   USE bibliography,                    ONLY: Golze2013,&
                                              Laino2005,&
                                              Laino2006,&
                                              cite_reference
   USE cell_methods,                    ONLY: write_cell
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_type,&
                                              get_cell
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_methods,               ONLY: create_small_subsys
   USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                              cp_subsys_type
   USE fist_environment,                ONLY: fist_init
   USE fist_environment_types,          ONLY: fist_env_create,&
                                              fist_env_get,&
                                              fist_env_set,&
                                              fist_environment_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
   USE global_types,                    ONLY: global_environment_type
   USE header,                          ONLY: qmmm_header
   USE input_constants,                 ONLY: &
        do_fist, do_multipole_section_off, do_multipole_section_on, do_qmmm, &
        do_qmmm_center_every_step, do_qmmm_center_never, do_qmmm_center_pbc_aware, &
        do_qmmm_center_setup_only, do_qmmm_none, do_qs
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE qmmm_init,                       ONLY: &
        assign_mm_charges_and_radius, move_or_add_atoms, print_image_charge_info, &
        print_qmmm_charges, print_qmmm_links, qmmm_init_gaussian_type, &
        qmmm_init_periodic_potential, qmmm_init_potential, setup_origin_mm_cell, setup_qmmm_links, &
        setup_qmmm_vars_mm, setup_qmmm_vars_qm
   USE qmmm_links_methods,              ONLY: qmmm_link_Imomm_coord
   USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
   USE qmmm_types,                      ONLY: qmmm_env_type
   USE qmmm_types_low,                  ONLY: add_set_release,&
                                              add_set_type,&
                                              add_shell_type,&
                                              qmmm_env_mm_create,&
                                              qmmm_env_mm_type,&
                                              qmmm_env_qm_create,&
                                              qmmm_env_qm_type,&
                                              qmmm_links_type
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_environment_type,&
                                              set_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_create'

   PUBLIC :: qmmm_env_create

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param force_env_section ...
!> \param qmmm_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param prev_subsys ...
!> \param ignore_outside_box ...
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qmmm_env_create(qmmm_env, root_section, para_env, globenv, &
                              force_env_section, qmmm_section, subsys_section, use_motion_section, prev_subsys, &
                              ignore_outside_box)
      TYPE(qmmm_env_type), INTENT(OUT)                   :: qmmm_env
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: force_env_section, qmmm_section, &
                                                            subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: prev_subsys
      LOGICAL, INTENT(in), OPTIONAL                      :: ignore_outside_box

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qmmm_env_create'

      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER                           :: qm_atom_type
      INTEGER                                            :: center_i, delta_charge, handle, iw, iw2, &
                                                            orig_charge, qmmm_coupl_type, &
                                                            use_multipole
      INTEGER, DIMENSION(:), POINTER                     :: mm_atom_index, mm_link_atoms, &
                                                            qm_atom_index
      LOGICAL                                            :: add_mm_charges, explicit, &
                                                            move_mm_charges, nocompatibility, &
                                                            qmmm_link, qmmm_link_Imomm, shell_model
      REAL(dp), DIMENSION(:), POINTER                    :: mm_atom_chrg, mm_el_pot_radius, &
                                                            mm_el_pot_radius_corr
      REAL(KIND=dp)                                      :: eps_mm_rspace
      REAL(KIND=dp), DIMENSION(3)                        :: abc_mm, abc_qm
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fist_scale_charge_link, &
                                                            mm_link_scale_factor
      TYPE(add_set_type), POINTER                        :: added_charges
      TYPE(add_shell_type), POINTER                      :: added_shells
      TYPE(cell_type), POINTER                           :: mm_cell, qm_cell_small, super_cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mm, subsys_qm
      TYPE(fist_environment_type), POINTER               :: fist_env
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qmmm_env_mm_type), POINTER                    :: qmmm_env_mm
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env_qm
      TYPE(qmmm_links_type), POINTER                     :: qmmm_links
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: multipole_section, print_gen, &
                                                            print_section, qmmm_periodic

      CALL timeset(routineN, handle)

      NULLIFY (qm_atom_index, mm_atom_index, qm_atom_type)
      NULLIFY (qmmm_env_qm, subsys_mm, subsys_qm, mm_cell, qm_cell_small)
      NULLIFY (mm_atom_chrg, mm_el_pot_radius, qmmm_env_mm, fist_env)
      NULLIFY (mm_link_atoms, mm_link_scale_factor, qmmm_links, added_charges, added_shells)
      NULLIFY (fist_scale_charge_link, print_section, fist_nonbond_env)
      NULLIFY (print_gen, logger, mm_el_pot_radius_corr, super_cell, pw_env)

      logger => cp_get_default_logger()

      ! citations
      CALL cite_reference(Laino2005)

      ! Input section...
      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF
      qmmm_periodic => section_vals_get_subs_vals(qmmm_section, "PERIODIC")
      multipole_section => section_vals_get_subs_vals(qmmm_section, "PERIODIC%MULTIPOLE")
      print_section => section_vals_get_subs_vals(qmmm_section, "PRINT")
      print_gen => section_vals_get_subs_vals(print_section, "PROGRAM_RUN_INFO")
      iw = cp_print_key_unit_nr(logger, print_gen, "", extension=".log")

      ! Create QM/MM Environments..
      ALLOCATE (qmmm_env_qm)
      CALL qmmm_env_qm_create(qmmm_env_qm)
      ALLOCATE (qmmm_env_mm)
      CALL qmmm_env_mm_create(qmmm_env_mm)

      ! Set up QM/MM Options
      CALL setup_qmmm_vars_mm(qmmm_section, &
                              qmmm_env_mm, &
                              qm_atom_index, &
                              mm_link_atoms, &
                              mm_link_scale_factor, &
                              fist_scale_charge_link, &
                              qmmm_coupl_type, &
                              qmmm_link)

      qmmm_env_mm%qm_atom_index => qm_atom_index
      qmmm_env_mm%mm_link_atoms => mm_link_atoms
      qmmm_env_mm%mm_link_scale_factor => mm_link_scale_factor
      qmmm_env_mm%fist_scale_charge_link => fist_scale_charge_link
      qmmm_env_mm%qmmm_coupl_type = qmmm_coupl_type
      qmmm_env_mm%qmmm_link = qmmm_link
      ! Center the qm subsys into the qm box
      CALL section_vals_val_get(qmmm_section, "CENTER", i_val=center_i)
      IF (center_i == do_qmmm_center_never) THEN
         qmmm_env_qm%center_qm_subsys = .FALSE.
         qmmm_env_qm%center_qm_subsys0 = .FALSE.
      ELSE IF (center_i == do_qmmm_center_setup_only) THEN
         qmmm_env_qm%center_qm_subsys = .FALSE.
         qmmm_env_qm%center_qm_subsys0 = .TRUE.
      ELSE IF (center_i == do_qmmm_center_every_step) THEN
         qmmm_env_qm%center_qm_subsys = .TRUE.
         qmmm_env_qm%center_qm_subsys0 = .TRUE.
      ELSE
         CPABORT("Unknown type of CENTER! ")
      END IF

      CALL section_vals_val_get(qmmm_section, "CENTER_TYPE", i_val=center_i)
      qmmm_env_qm%center_qm_subsys_pbc_aware = (center_i == do_qmmm_center_pbc_aware)

      ! Compatibility with the QM/MM in CPMD code
      CALL section_vals_val_get(qmmm_section, "NOCOMPATIBILITY", l_val=nocompatibility)
      qmmm_env_qm%compatibility = .NOT. nocompatibility

      ! Parallel scheme for the long range
      CALL section_vals_val_get(qmmm_section, "PARALLEL_SCHEME", &
                                i_val=qmmm_env_qm%par_scheme)

      ! Periodic boundary condition calculation
      CALL section_vals_get(qmmm_periodic, explicit=explicit)
      qmmm_env_qm%periodic = explicit
      !multipole section is switched on by default; switched off only if explicitly stated
      IF (qmmm_env_qm%periodic) qmmm_env_qm%multipole = .TRUE.
      CALL section_vals_get(multipole_section, explicit=explicit)
      CALL section_vals_val_get(multipole_section, "_SECTION_PARAMETERS_", i_val=use_multipole)
      IF (explicit .AND. use_multipole == do_multipole_section_off) qmmm_env_qm%multipole = .FALSE.
      IF (explicit .AND. use_multipole == do_multipole_section_on) qmmm_env_qm%multipole = .TRUE.
      IF (qmmm_env_qm%periodic .AND. qmmm_env_qm%multipole) CALL cite_reference(Laino2006)
      IF (qmmm_coupl_type == do_qmmm_none) THEN
         IF (qmmm_env_qm%periodic) &
            CALL cp_warn(__LOCATION__, &
                         "QMMM periodic calculation with coupling NONE was requested! "// &
                         "Switching off the periodic keyword since periodic and non-periodic "// &
                         "calculation with coupling NONE represent the same method! ")
         qmmm_env_qm%periodic = .FALSE.
      END IF

      ! First Initialize Fist...
      CALL section_vals_val_set(force_env_section, "METHOD", i_val=do_fist)
      ALLOCATE (fist_env)
      CALL fist_env_create(fist_env, para_env=para_env)
      CALL fist_env_set(fist_env, qmmm=.TRUE., qmmm_env=qmmm_env_mm)
      CALL fist_init(fist_env, root_section, para_env, force_env_section, &
                     subsys_section, use_motion_section, prev_subsys=prev_subsys)
      CALL fist_env_get(fist_env, subsys=subsys_mm, cell=mm_cell)
      mm_cell%tag = "CELL_MM"

      ! Set up QM/MM Options
      CALL setup_qmmm_vars_qm(qmmm_section, &
                              qmmm_env_qm, &
                              subsys_mm, &
                              qm_atom_type, &
                              qm_atom_index, &
                              mm_atom_index, &
                              qm_cell_small, &
                              qmmm_coupl_type, &
                              eps_mm_rspace, &
                              qmmm_link, &
                              para_env)

      qmmm_env_qm%qm_atom_index => qm_atom_index
      qmmm_env_qm%mm_atom_index => mm_atom_index
      qmmm_env_qm%eps_mm_rspace = eps_mm_rspace
      qmmm_env_qm%qmmm_coupl_type = qmmm_coupl_type
      qmmm_env_qm%qmmm_link = qmmm_link
      qmmm_env_qm%num_qm_atoms = SIZE(qm_atom_index)
      qmmm_env_qm%num_mm_atoms = SIZE(mm_atom_index)
      IF (qmmm_env_qm%image_charge) THEN
         qmmm_env_qm%num_image_mm_atoms = SIZE(qmmm_env_qm%image_charge_pot%image_mm_list)
         CALL cite_reference(Golze2013)
      END IF

      ! Duplicate structure for link atoms
      IF (qmmm_link) THEN
         IF (ASSOCIATED(mm_link_atoms)) THEN
            ALLOCATE (qmmm_env_qm%mm_link_atoms(SIZE(mm_link_atoms)))
            qmmm_env_qm%mm_link_atoms = mm_link_atoms
         END IF
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(A,I26)') " Number of QM atoms: ", qmmm_env_qm%num_qm_atoms
         WRITE (iw, '(A,I26)') " Number of MM atoms: ", qmmm_env_qm%num_mm_atoms
         IF (qmmm_env_qm%image_charge) THEN
            WRITE (iw, '(A,I8)') " Number of MM atoms with image charge: ", &
               qmmm_env_qm%num_image_mm_atoms
         END IF
         CALL write_cell(mm_cell, subsys_section)
      END IF
      CALL get_cell(qm_cell_small, abc=abc_qm)
      CALL get_cell(mm_cell, abc=abc_mm)

      IF (qmmm_env_qm%image_charge) THEN
         IF (ANY(ABS(abc_mm - abc_qm) > 1.0E-12)) &
            CPABORT("QM and MM box need to have the same size when using image charges")
      END IF

      ! Assign charges and mm_el_pot_radius from fist_topology
      CALL fist_env_get(fist_env, fist_nonbond_env=fist_nonbond_env)
      ALLOCATE (mm_atom_chrg(SIZE(mm_atom_index)))
      ALLOCATE (mm_el_pot_radius(SIZE(mm_atom_index)))
      ALLOCATE (mm_el_pot_radius_corr(SIZE(mm_atom_index)))
      mm_atom_chrg = 0.0_dp
      mm_el_pot_radius = 0.0_dp
      mm_el_pot_radius_corr = 0.0_dp

      CALL assign_mm_charges_and_radius(subsys=subsys_mm, &
                                        charges=fist_nonbond_env%charges, &
                                        mm_atom_chrg=mm_atom_chrg, &
                                        mm_el_pot_radius=mm_el_pot_radius, &
                                        mm_el_pot_radius_corr=mm_el_pot_radius_corr, &
                                        mm_atom_index=mm_atom_index, &
                                        mm_link_atoms=mm_link_atoms, &
                                        mm_link_scale_factor=mm_link_scale_factor, &
                                        added_shells=added_shells, &
                                        shell_model=shell_model)

      qmmm_env_qm%mm_atom_chrg => mm_atom_chrg
      qmmm_env_qm%mm_el_pot_radius => mm_el_pot_radius
      qmmm_env_qm%mm_el_pot_radius_corr => mm_el_pot_radius_corr
      qmmm_env_qm%added_shells => added_shells

      qmmm_link_Imomm = .FALSE.
      IF (qmmm_link) THEN
         CALL setup_qmmm_links(qmmm_section, qmmm_links, mm_el_pot_radius, &
                               mm_el_pot_radius_corr, mm_atom_index, iw)
         qmmm_env_qm%qmmm_links => qmmm_links

         CALL print_qmmm_links(qmmm_section, qmmm_links)

         CALL add_set_release(qmmm_env_qm%added_charges)
         CALL move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                                mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr, &
                                added_charges, mm_atom_index)
         qmmm_env_qm%move_mm_charges = move_mm_charges
         qmmm_env_qm%add_mm_charges = add_mm_charges
         qmmm_env_qm%added_charges => added_charges
         IF (ASSOCIATED(qmmm_links%imomm)) qmmm_link_imomm = (SIZE(qmmm_links%imomm) /= 0)
      END IF

      CALL print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius, &
                              mm_el_pot_radius_corr, qmmm_env_qm%added_charges, &
                              qmmm_env_qm%added_shells, qmmm_section, nocompatibility, shell_model)
      IF (qmmm_env_qm%image_charge) THEN
         CALL print_image_charge_info(qmmm_env_qm, qmmm_section)
      END IF

      CALL section_vals_val_get(qmmm_section, "DELTA_CHARGE", i_val=delta_charge)
      CALL section_vals_val_get(force_env_section, "DFT%CHARGE", i_val=orig_charge)
      CALL section_vals_val_set(force_env_section, "DFT%CHARGE", i_val=orig_charge + delta_charge)

      CALL section_vals_val_set(force_env_section, "METHOD", i_val=do_qs)
      CALL create_small_subsys(subsys_qm, &
                               big_subsys=subsys_mm, small_para_env=para_env, &
                               small_cell=qm_cell_small, sub_atom_index=qm_atom_index, &
                               sub_atom_kind_name=qm_atom_type, para_env=para_env, &
                               force_env_section=force_env_section, subsys_section=subsys_section, &
                               ignore_outside_box=ignore_outside_box)
      IF (qmmm_link_imomm) CALL qmmm_link_Imomm_coord(qmmm_links, subsys_qm%particles%els, &
                                                      qm_atom_index)
      ALLOCATE (qs_env)
      CALL qs_env_create(qs_env, globenv)
      CALL qs_init(qs_env, para_env, root_section, globenv=globenv, cp_subsys=subsys_qm, &
                   cell=qm_cell_small, qmmm=.TRUE., qmmm_env_qm=qmmm_env_qm, &
                   force_env_section=force_env_section, subsys_section=subsys_section, &
                   use_motion_section=use_motion_section)
      CALL cp_subsys_release(subsys_qm)

      IF (qmmm_env_qm%periodic) THEN
         IF (.NOT. ASSOCIATED(super_cell)) THEN
            ALLOCATE (super_cell)
         END IF
         CALL cell_clone(mm_cell, super_cell, tag="SUPER_CELL")
         CALL set_qs_env(qs_env, super_cell=super_cell, qmmm_periodic=qmmm_env_qm%periodic)
         CALL cell_release(super_cell)
      END IF
      CALL section_vals_val_set(force_env_section, "DFT%CHARGE", i_val=orig_charge)
      CALL cp_print_key_finished_output(iw, logger, print_gen, "")
      iw2 = cp_print_key_unit_nr(logger, qmmm_section, "PRINT%PROGRAM_BANNER", &
                                 extension=".qmmmLog")
      CALL qmmm_header(iw2)
      CALL cp_print_key_finished_output(iw2, logger, qmmm_section, &
                                        "PRINT%PROGRAM_BANNER")
      !
      ! Initialize MM Potential fitted with Gaussian
      !
      CALL qmmm_init_gaussian_type(qmmm_env_qm=qmmm_env_qm, &
                                   para_env=para_env, &
                                   qs_env=qs_env, &
                                   mm_atom_chrg=mm_atom_chrg, &
                                   added_charges=qmmm_env_qm%added_charges, &
                                   added_shells=qmmm_env_qm%added_shells, &
                                   print_section=print_section, &
                                   qmmm_section=qmmm_section)
      !
      ! Initialize the MM potential stored on vector
      !
      CALL qmmm_init_potential(qmmm_env_qm=qmmm_env_qm, &
                               mm_cell=mm_cell, &
                               added_charges=qmmm_env_qm%added_charges, &
                               added_shells=qmmm_env_qm%added_shells, &
                               print_section=print_section)
      !
      ! Initialize the qmmm_pw_grid
      !
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL qmmm_pw_grid_init(qmmm_env=qmmm_env_qm, &
                             pw_env=pw_env)
      !
      ! Initialize the MM periodic potential
      !
      CALL qmmm_init_periodic_potential(qmmm_env_qm=qmmm_env_qm, &
                                        qm_cell_small=qm_cell_small, &
                                        mm_cell=mm_cell, &
                                        para_env=para_env, &
                                        qs_env=qs_env, &
                                        added_charges=qmmm_env_qm%added_charges, &
                                        added_shells=qmmm_env_qm%added_shells, &
                                        qmmm_periodic=qmmm_periodic, &
                                        print_section=print_section, &
                                        mm_atom_chrg=mm_atom_chrg)
      !
      ! Preparing for PBC...
      !
      CALL setup_origin_mm_cell(qmmm_section, qmmm_env_qm, qm_cell_small, &
                                dr=pw_env%pw_pools(pw_env%auxbas_grid)%pool%pw_grid%dr)

      CALL cell_release(qm_cell_small)

      ! assemble the actual qmmm_env
      qmmm_env%qs_env => qs_env
      qmmm_env%fist_env => fist_env
      qmmm_env%qm => qmmm_env_qm

      CALL section_vals_val_set(force_env_section, "METHOD", i_val=do_qmmm)
      DEALLOCATE (qm_atom_type)

      CALL timestop(handle)

   END SUBROUTINE qmmm_env_create

END MODULE qmmm_create
