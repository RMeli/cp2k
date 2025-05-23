!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Handles the MPI communication of the swarm framework.
!> \author Ole Schuett
! **************************************************************************************************
MODULE swarm_mpi
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_iter_types,                   ONLY: cp_iteration_info_create,&
                                              cp_iteration_info_release,&
                                              cp_iteration_info_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_any_source,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE swarm_message,                   ONLY: swarm_message_get,&
                                              swarm_message_mpi_bcast,&
                                              swarm_message_mpi_recv,&
                                              swarm_message_mpi_send,&
                                              swarm_message_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_mpi'

   PUBLIC :: swarm_mpi_type, swarm_mpi_init, swarm_mpi_finalize
   PUBLIC :: swarm_mpi_send_report, swarm_mpi_recv_report
   PUBLIC :: swarm_mpi_send_command, swarm_mpi_recv_command

   TYPE swarm_mpi_type
      TYPE(mp_para_env_type), POINTER          :: world => Null()
      TYPE(mp_para_env_type), POINTER          :: worker => Null()
      TYPE(mp_para_env_type), POINTER          :: master => Null()
      INTEGER, DIMENSION(:), ALLOCATABLE       :: wid2group
      CHARACTER(LEN=default_path_length)       :: master_output_path = ""
   END TYPE swarm_mpi_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize MPI communicators for a swarm run.
!> \param swarm_mpi ...
!> \param world_para_env ...
!> \param root_section ...
!> \param n_workers ...
!> \param worker_id ...
!> \param iw ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_init(swarm_mpi, world_para_env, root_section, n_workers, worker_id, iw)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(mp_para_env_type), POINTER                    :: world_para_env
      TYPE(section_vals_type), POINTER                   :: root_section
      INTEGER, INTENT(IN)                                :: n_workers
      INTEGER, INTENT(OUT)                               :: worker_id
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: n_groups_created, pe_per_worker, &
                                                            subgroup_rank, subgroup_size
      TYPE(mp_comm_type)                                 :: subgroup
      LOGICAL                                            :: im_the_master
      INTEGER, DIMENSION(:), POINTER                     :: group_distribution_p
      INTEGER, DIMENSION(0:world_para_env%num_pe-2), &
         TARGET                                          :: group_distribution

! ====== Setup of MPI-Groups ======

      worker_id = -1
      swarm_mpi%world => world_para_env

      IF (MOD(swarm_mpi%world%num_pe - 1, n_workers) /= 0) THEN
         CPABORT("number of processors-1 is not divisible by n_workers.")
      END IF
      IF (swarm_mpi%world%num_pe < n_workers + 1) THEN
         CPABORT("There are not enough processes for n_workers + 1. Aborting.")
      END IF

      pe_per_worker = (swarm_mpi%world%num_pe - 1)/n_workers

      IF (iw > 0) THEN
         WRITE (iw, '(A,45X,I8)') " SWARM| Number of mpi ranks", swarm_mpi%world%num_pe
         WRITE (iw, '(A,47X,I8)') " SWARM| Number of workers", n_workers
      END IF

      ! the last task becomes the master. Preseves node-alignment of other tasks.
      im_the_master = (swarm_mpi%world%mepos == swarm_mpi%world%num_pe - 1)

      ! First split split para_env into a master- and a workers-groups...
      IF (im_the_master) THEN
         ALLOCATE (swarm_mpi%master)
         CALL swarm_mpi%master%from_split(swarm_mpi%world, 1)
         subgroup_size = swarm_mpi%master%num_pe
         subgroup_rank = swarm_mpi%master%mepos
         IF (swarm_mpi%master%num_pe /= 1) CPABORT("mp_comm_split_direct failed (master)")
      ELSE
         CALL subgroup%from_split(swarm_mpi%world, 2)
         subgroup_size = subgroup%num_pe
         subgroup_rank = subgroup%mepos
         IF (subgroup_size /= swarm_mpi%world%num_pe - 1) CPABORT("mp_comm_split_direct failed (worker)")
      END IF

      ALLOCATE (swarm_mpi%wid2group(n_workers))
      swarm_mpi%wid2group = 0

      IF (.NOT. im_the_master) THEN
         ! ...then split workers-group into n_workers groups - one for each worker.
         group_distribution_p => group_distribution
         ALLOCATE (swarm_mpi%worker)
         CALL swarm_mpi%worker%from_split(subgroup, n_groups_created, group_distribution_p, n_subgroups=n_workers)
         worker_id = group_distribution(subgroup_rank) + 1 ! shall start by 1
         IF (n_groups_created /= n_workers) CPABORT("mp_comm_split failed.")
         CALL subgroup%free()

         !WRITE (*,*) "this is worker ", worker_id, swarm_mpi%worker%mepos, swarm_mpi%worker%num_pe

         ! collect world-ranks of each worker groups rank-0 node
         IF (swarm_mpi%worker%mepos == 0) &
            swarm_mpi%wid2group(worker_id) = swarm_mpi%world%mepos

      END IF

      CALL swarm_mpi%world%sum(swarm_mpi%wid2group)
      !WRITE (*,*), "wid2group table: ",swarm_mpi%wid2group

      CALL logger_init_master(swarm_mpi)
      CALL logger_init_worker(swarm_mpi, root_section, worker_id)
   END SUBROUTINE swarm_mpi_init

! **************************************************************************************************
!> \brief Helper routine for swarm_mpi_init, configures the master's logger.
!> \param swarm_mpi ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE logger_init_master(swarm_mpi)
      TYPE(swarm_mpi_type)                               :: swarm_mpi

      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger

! broadcast master_output_path to all ranks

      IF (swarm_mpi%world%is_source()) THEN
         logger => cp_get_default_logger()
         output_unit = logger%default_local_unit_nr
         swarm_mpi%master_output_path = output_unit2path(output_unit)
         IF (output_unit /= default_output_unit) &
            CLOSE (output_unit)
      END IF

      CALL swarm_mpi%world%bcast(swarm_mpi%master_output_path)

      IF (ASSOCIATED(swarm_mpi%master)) &
         CALL error_add_new_logger(swarm_mpi%master, swarm_mpi%master_output_path)
   END SUBROUTINE logger_init_master

! **************************************************************************************************
!> \brief Helper routine for logger_init_master, inquires filename for given unit.
!> \param output_unit ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION output_unit2path(output_unit) RESULT(output_path)
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=default_path_length)                 :: output_path

      output_path = "__STD_OUT__"
      IF (output_unit /= default_output_unit) &
         INQUIRE (unit=output_unit, name=output_path)
   END FUNCTION output_unit2path

! **************************************************************************************************
!> \brief Helper routine for swarm_mpi_init, configures the workers's logger.
!> \param swarm_mpi ...
!> \param root_section ...
!> \param worker_id ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE logger_init_worker(swarm_mpi, root_section, worker_id)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(section_vals_type), POINTER                   :: root_section
      INTEGER                                            :: worker_id

      CHARACTER(LEN=default_path_length)                 :: output_path
      CHARACTER(len=default_string_length)               :: new_project_name, project_name, &
                                                            worker_name
      TYPE(cp_iteration_info_type), POINTER              :: new_iter_info
      TYPE(cp_logger_type), POINTER                      :: old_logger

      NULLIFY (old_logger, new_iter_info)
      IF (ASSOCIATED(swarm_mpi%worker)) THEN
         old_logger => cp_get_default_logger()
         project_name = old_logger%iter_info%project_name
         IF (worker_id > 99999) THEN
            CPABORT("Did not expect so many workers.")
         END IF
         WRITE (worker_name, "(A,I5.5)") 'WORKER', worker_id
         IF (LEN_TRIM(project_name) + 1 + LEN_TRIM(worker_name) > default_string_length) THEN
            CPABORT("project name too long")
         END IF
         output_path = TRIM(project_name)//"-"//TRIM(worker_name)//".out"
         new_project_name = TRIM(project_name)//"-"//TRIM(worker_name)
         CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name)
         CALL cp_iteration_info_create(new_iter_info, new_project_name)
         CALL error_add_new_logger(swarm_mpi%worker, output_path, new_iter_info)
         CALL cp_iteration_info_release(new_iter_info)
      END IF
   END SUBROUTINE logger_init_worker

! **************************************************************************************************
!> \brief Helper routine for logger_init_master and logger_init_worker
!> \param para_env ...
!> \param output_path ...
!> \param iter_info ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE error_add_new_logger(para_env, output_path, iter_info)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=default_path_length)                 :: output_path
      TYPE(cp_iteration_info_type), OPTIONAL, POINTER    :: iter_info

      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: new_logger, old_logger

      NULLIFY (new_logger, old_logger)
      output_unit = -1
      IF (para_env%is_source()) THEN
         ! open output_unit according to output_path
         output_unit = default_output_unit
         IF (output_path /= "__STD_OUT__") &
            CALL open_file(file_name=output_path, file_status="UNKNOWN", &
                           file_action="WRITE", file_position="APPEND", unit_number=output_unit)
      END IF

      old_logger => cp_get_default_logger()
      CALL cp_logger_create(new_logger, para_env=para_env, &
                            default_global_unit_nr=output_unit, close_global_unit_on_dealloc=.FALSE., &
                            template_logger=old_logger, iter_info=iter_info)

      CALL cp_add_default_logger(new_logger)
      CALL cp_logger_release(new_logger)
   END SUBROUTINE error_add_new_logger

! **************************************************************************************************
!> \brief Finalizes the MPI communicators of a swarm run.
!> \param swarm_mpi ...
!> \param root_section ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_finalize(swarm_mpi, root_section)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(section_vals_type), POINTER                   :: root_section

      CALL swarm_mpi%world%sync()
      CALL logger_finalize(swarm_mpi, root_section)

      IF (ASSOCIATED(swarm_mpi%worker)) CALL mp_para_env_release(swarm_mpi%worker)
      IF (ASSOCIATED(swarm_mpi%master)) CALL mp_para_env_release(swarm_mpi%master)
      NULLIFY (swarm_mpi%worker, swarm_mpi%master)
      DEALLOCATE (swarm_mpi%wid2group)
   END SUBROUTINE swarm_mpi_finalize

! **************************************************************************************************
!> \brief Helper routine for swarm_mpi_finalize, restores the original loggers
!> \param swarm_mpi ...
!> \param root_section ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE logger_finalize(swarm_mpi, root_section)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(section_vals_type), POINTER                   :: root_section

      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger, old_logger

      NULLIFY (logger, old_logger)
      logger => cp_get_default_logger()
      output_unit = logger%default_local_unit_nr
      IF (output_unit > 0 .AND. output_unit /= default_output_unit) &
         CALL close_file(output_unit)

      CALL cp_rm_default_logger() !pops the top-most logger
      old_logger => cp_get_default_logger()

      ! restore GLOBAL%PROJECT_NAME
      CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", &
                                c_val=old_logger%iter_info%project_name)

      CALL swarm_mpi%world%sync()

      ! do this only on master's rank 0
      IF (swarm_mpi%world%is_source() .AND. output_unit /= default_output_unit) THEN
         output_unit = old_logger%default_local_unit_nr
         OPEN (unit=output_unit, file=swarm_mpi%master_output_path, &
               status="UNKNOWN", action="WRITE", position="APPEND")
      END IF
   END SUBROUTINE logger_finalize

! **************************************************************************************************
!> \brief Sends a report via MPI
!> \param swarm_mpi ...
!> \param report ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_send_report(swarm_mpi, report)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(swarm_message_type)                           :: report

      INTEGER                                            :: dest, tag

! Only rank-0 of worker group sends its report

      IF (swarm_mpi%worker%is_source()) THEN
         dest = swarm_mpi%world%num_pe - 1
         tag = 42
         CALL swarm_message_mpi_send(report, group=swarm_mpi%world, dest=dest, tag=tag)
      END IF

   END SUBROUTINE swarm_mpi_send_report

! **************************************************************************************************
!> \brief Receives a report via MPI
!> \param swarm_mpi ...
!> \param report ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_recv_report(swarm_mpi, report)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(swarm_message_type), INTENT(OUT)              :: report

      INTEGER                                            :: src, tag

      tag = 42
      src = mp_any_source

      CALL swarm_message_mpi_recv(report, group=swarm_mpi%world, src=src, tag=tag)

   END SUBROUTINE swarm_mpi_recv_report

! **************************************************************************************************
!> \brief Sends a command via MPI
!> \param swarm_mpi ...
!> \param cmd ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_send_command(swarm_mpi, cmd)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(swarm_message_type)                           :: cmd

      INTEGER                                            :: dest, tag, worker_id

      CALL swarm_message_get(cmd, "worker_id", worker_id)
      tag = 42
      dest = swarm_mpi%wid2group(worker_id)

      CALL swarm_message_mpi_send(cmd, group=swarm_mpi%world, dest=dest, tag=tag)

   END SUBROUTINE swarm_mpi_send_command

! **************************************************************************************************
!> \brief Receives a command via MPI and broadcasts it within a worker.
!> \param swarm_mpi ...
!> \param cmd ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_mpi_recv_command(swarm_mpi, cmd)
      TYPE(swarm_mpi_type)                               :: swarm_mpi
      TYPE(swarm_message_type), INTENT(OUT)              :: cmd

      INTEGER                                            :: src, tag

! This is a two step communication schema.
! First: The rank-0 of the worker groups receives the command from the master.

      IF (swarm_mpi%worker%is_source()) THEN
         src = swarm_mpi%world%num_pe - 1 !
         tag = 42
         CALL swarm_message_mpi_recv(cmd, group=swarm_mpi%world, src=src, tag=tag)

      END IF

!     ! Second: The command is broadcasted within the worker group.
      CALL swarm_message_mpi_bcast(cmd, src=swarm_mpi%worker%source, group=swarm_mpi%worker)

   END SUBROUTINE swarm_mpi_recv_command

END MODULE swarm_mpi

