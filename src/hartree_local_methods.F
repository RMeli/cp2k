!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE hartree_local_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE hartree_local_types,             ONLY: allocate_ecoul_1center,&
                                              ecoul_1center_type,&
                                              hartree_local_type,&
                                              set_ecoul_1c
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_periodic,&
                                              pw_poisson_type
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_local_rho_types,              ONLY: get_local_rho,&
                                              local_rho_type,&
                                              rhoz_type
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              rho0_atom_type,&
                                              rho0_mpole_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hartree_local_methods'

   ! Public Subroutine

   PUBLIC :: init_coulomb_local, calculate_Vh_1center, Vh_1c_gg_integrals

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param hartree_local ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE init_coulomb_local(hartree_local, natom)

      TYPE(hartree_local_type), POINTER                  :: hartree_local
      INTEGER, INTENT(IN)                                :: natom

      CHARACTER(len=*), PARAMETER :: routineN = 'init_coulomb_local'

      INTEGER                                            :: handle
      TYPE(ecoul_1center_type), DIMENSION(:), POINTER    :: ecoul_1c

      CALL timeset(routineN, handle)

      NULLIFY (ecoul_1c)
      !   Allocate and Initialize 1-center Potentials and Integrals
      CALL allocate_ecoul_1center(ecoul_1c, natom)
      hartree_local%ecoul_1c => ecoul_1c

      CALL timestop(handle)

   END SUBROUTINE init_coulomb_local

! **************************************************************************************************
!> \brief Calculates Hartree potential for hard and soft densities (including
!>        nuclear charge and compensation charges) using numerical integration
!> \param vrad_h ...
!> \param vrad_s ...
!> \param rrad_h ...
!> \param rrad_s ...
!> \param rrad_0 ...
!> \param rrad_z ...
!> \param grid_atom ...
!> \par History
!>      05.2012 JGH refactoring
!> \author ??
! **************************************************************************************************
   SUBROUTINE calculate_Vh_1center(vrad_h, vrad_s, rrad_h, rrad_s, rrad_0, rrad_z, grid_atom)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: vrad_h, vrad_s
      TYPE(rho_atom_coeff), DIMENSION(:), INTENT(IN)     :: rrad_h, rrad_s
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rrad_0
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: rrad_z
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_Vh_1center'

      INTEGER                                            :: handle, ir, iso, ispin, l_ang, &
                                                            max_s_harm, nchannels, nr, nspins
      REAL(dp)                                           :: I1_down, I1_up, I2_down, I2_up, prefactor
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rho_1, rho_2
      REAL(dp), DIMENSION(:), POINTER                    :: wr
      REAL(dp), DIMENSION(:, :), POINTER                 :: oor2l, r2l

      CALL timeset(routineN, handle)

      nr = grid_atom%nr
      max_s_harm = SIZE(vrad_h, 2)
      nspins = SIZE(rrad_h, 1)
      nchannels = SIZE(rrad_0, 2)

      r2l => grid_atom%rad2l
      oor2l => grid_atom%oorad2l
      wr => grid_atom%wr

      ALLOCATE (rho_1(nr, max_s_harm), rho_2(nr, max_s_harm))
      rho_1 = 0.0_dp
      rho_2 = 0.0_dp

      !   Case lm = 0
      rho_1(:, 1) = rrad_z(:)
      rho_2(:, 1) = rrad_0(:, 1)

      DO iso = 2, nchannels
         rho_2(:, iso) = rrad_0(:, iso)
      END DO

      DO iso = 1, max_s_harm
         DO ispin = 1, nspins
            rho_1(:, iso) = rho_1(:, iso) + rrad_h(ispin)%r_coef(:, iso)
            rho_2(:, iso) = rho_2(:, iso) + rrad_s(ispin)%r_coef(:, iso)
         END DO

         l_ang = indso(1, iso)
         prefactor = fourpi/(2._dp*l_ang + 1._dp)

         rho_1(:, iso) = rho_1(:, iso)*wr(:)
         rho_2(:, iso) = rho_2(:, iso)*wr(:)

         I1_up = 0.0_dp
         I1_down = 0.0_dp
         I2_up = 0.0_dp
         I2_down = 0.0_dp

         I1_up = r2l(nr, l_ang)*rho_1(nr, iso)
         I2_up = r2l(nr, l_ang)*rho_2(nr, iso)

         DO ir = nr - 1, 1, -1
            I1_down = I1_down + oor2l(ir, l_ang + 1)*rho_1(ir, iso)
            I2_down = I2_down + oor2l(ir, l_ang + 1)*rho_2(ir, iso)
         END DO

         vrad_h(nr, iso) = vrad_h(nr, iso) + prefactor* &
                           (oor2l(nr, l_ang + 1)*I1_up + r2l(nr, l_ang)*I1_down)
         vrad_s(nr, iso) = vrad_s(nr, iso) + prefactor* &
                           (oor2l(nr, l_ang + 1)*I2_up + r2l(nr, l_ang)*I2_down)

         DO ir = nr - 1, 1, -1
            I1_up = I1_up + r2l(ir, l_ang)*rho_1(ir, iso)
            I1_down = I1_down - oor2l(ir, l_ang + 1)*rho_1(ir, iso)
            I2_up = I2_up + r2l(ir, l_ang)*rho_2(ir, iso)
            I2_down = I2_down - oor2l(ir, l_ang + 1)*rho_2(ir, iso)

            vrad_h(ir, iso) = vrad_h(ir, iso) + prefactor* &
                              (oor2l(ir, l_ang + 1)*I1_up + r2l(ir, l_ang)*I1_down)
            vrad_s(ir, iso) = vrad_s(ir, iso) + prefactor* &
                              (oor2l(ir, l_ang + 1)*I2_up + r2l(ir, l_ang)*I2_down)

         END DO

      END DO

      DEALLOCATE (rho_1, rho_2)

      CALL timestop(handle)

   END SUBROUTINE calculate_Vh_1center

! **************************************************************************************************
!> \brief Calculates one center GAPW Hartree energies and matrix elements
!>        Hartree potentials are input
!>        Takes possible background charge into account
!>        Special case for densities without core charge
!> \param qs_env ...
!> \param energy_hartree_1c ...
!> \param ecoul_1c ...
!> \param local_rho_set ...
!> \param para_env ...
!> \param tddft ...
!> \param local_rho_set_2nd ...
!> \param core_2nd ...
!> \par History
!>      05.2012 JGH refactoring
!> \author ??
! **************************************************************************************************
   SUBROUTINE Vh_1c_gg_integrals(qs_env, energy_hartree_1c, ecoul_1c, local_rho_set, para_env, tddft, local_rho_set_2nd, &
                                 core_2nd)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(kind=dp), INTENT(out)                         :: energy_hartree_1c
      TYPE(ecoul_1center_type), DIMENSION(:), POINTER    :: ecoul_1c
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: tddft
      TYPE(local_rho_type), OPTIONAL, POINTER            :: local_rho_set_2nd
      LOGICAL, INTENT(IN), OPTIONAL                      :: core_2nd

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Vh_1c_gg_integrals'

      INTEGER :: bo(2), handle, iat, iatom, ikind, ipgf1, is1, iset1, iso, l_ang, llmax, lmax0, &
         lmax0_2nd, lmax_0, m1, max_iso, max_iso_not0, max_s_harm, maxl, maxso, mepos, n1, nat, &
         nchan_0, nkind, nr, nset, nsotot, nspins, num_pe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf
      LOGICAL                                            :: core_charge, l_2nd_local_rho, &
                                                            my_core_2nd, my_periodic, paw_atom
      REAL(dp)                                           :: back_ch, factor
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gexp, sqrtwr
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aVh1b_00, aVh1b_hh, aVh1b_ss, g0_h_w
      REAL(dp), DIMENSION(:), POINTER                    :: rrad_z, vrrad_z
      REAL(dp), DIMENSION(:, :), POINTER                 :: g0_h, g0_h_2nd, gsph, rrad_0, Vh1_h, &
                                                            Vh1_s, vrrad_0, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG, Qlm_gg, Qlm_gg_2nd
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rho0_atom_type), DIMENSION(:), POINTER        :: rho0_atom_set, rho0_atom_set_2nd
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole, rho0_mpole_2nd
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set, rho_atom_set_2nd
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(rhoz_type), DIMENSION(:), POINTER             :: rhoz_set, rhoz_set_2nd

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, poisson_env, pw_env, qs_charges)
      NULLIFY (atomic_kind_set, qs_kind_set, rho_atom_set, rho0_atom_set)
      NULLIFY (rho0_mpole, rhoz_set)
      NULLIFY (atom_list, grid_atom, harmonics)
      NULLIFY (basis_1c, lmin, lmax, npgf, zet)
      NULLIFY (gsph)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      pw_env=pw_env, qs_charges=qs_charges)

      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      my_periodic = (poisson_env%method == pw_poisson_periodic)

      back_ch = qs_charges%background*cell%deth

      ! rhoz_set is not accessed in TDDFT
      CALL get_local_rho(local_rho_set, rho_atom_set, rho0_atom_set, rho0_mpole, rhoz_set) ! for integral space

      ! for forces we need a second local_rho_set
      l_2nd_local_rho = .FALSE.
      IF (PRESENT(local_rho_set_2nd)) THEN ! for potential
         l_2nd_local_rho = .TRUE.
         NULLIFY (rho_atom_set_2nd, rho0_atom_set_2nd, rhoz_set_2nd) ! for potential
         CALL get_local_rho(local_rho_set_2nd, rho_atom_set_2nd, rho0_atom_set_2nd, rho0_mpole_2nd, rhoz_set=rhoz_set_2nd)
      END IF

      nkind = SIZE(atomic_kind_set, 1)
      nspins = dft_control%nspins

      core_charge = .NOT. tddft  ! for forces mixed version
      my_core_2nd = .TRUE.
      IF (PRESENT(core_2nd)) my_core_2nd = .NOT. core_2nd   ! if my_core_2nd true, include core charge

      ! The aim of the following code was to return immediately if the subroutine
      ! was called for triplet excited states in spin-restricted case. This check
      ! is also performed before invocation of this subroutine. It should be save
      ! to remove the optional argument 'do_triplet' from the subroutine interface.
      !IF (tddft) THEN
      !   CPASSERT(PRESENT(do_triplet))
      !   IF (nspins == 1 .AND. do_triplet) RETURN
      !END IF

      CALL get_qs_kind_set(qs_kind_set, maxg_iso_not0=max_iso)
      CALL get_rho0_mpole(rho0_mpole=rho0_mpole, lmax_0=lmax_0)

      !   Put to 0 the local hartree energy contribution from 1 center integrals
      energy_hartree_1c = 0.0_dp

      !   Here starts the loop over all the atoms
      DO ikind = 1, nkind

         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          grid_atom=grid_atom, &
                          harmonics=harmonics, ngrid_rad=nr, &
                          max_iso_not0=max_iso_not0, paw_atom=paw_atom)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=basis_1c, basis_type="GAPW_1C")

         IF (paw_atom) THEN
!===========    PAW   ===============
            CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                                   maxso=maxso, npgf=npgf, maxl=maxl, &
                                   nset=nset, zet=zet)

            max_s_harm = harmonics%max_s_harm
            llmax = harmonics%llmax

            nsotot = maxso*nset
            ALLOCATE (gsph(nr, nsotot))
            ALLOCATE (gexp(nr))
            ALLOCATE (sqrtwr(nr), g0_h_w(nr, 0:lmax_0))

            NULLIFY (Vh1_h, Vh1_s)
            ALLOCATE (Vh1_h(nr, max_iso_not0))
            ALLOCATE (Vh1_s(nr, max_iso_not0))

            ALLOCATE (aVh1b_hh(nsotot, nsotot))
            ALLOCATE (aVh1b_ss(nsotot, nsotot))
            ALLOCATE (aVh1b_00(nsotot, nsotot))
            ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

            NULLIFY (Qlm_gg, g0_h)
            CALL get_rho0_mpole(rho0_mpole=rho0_mpole, ikind=ikind, &
                                l0_ikind=lmax0, &
                                Qlm_gg=Qlm_gg, g0_h=g0_h) ! Qlm_gg of density

            IF (PRESENT(local_rho_set_2nd)) THEN ! for potential
               NULLIFY (Qlm_gg_2nd, g0_h_2nd)
               CALL get_rho0_mpole(rho0_mpole=rho0_mpole_2nd, ikind=ikind, &
                                   l0_ikind=lmax0_2nd, &
                                   Qlm_gg=Qlm_gg_2nd, g0_h=g0_h_2nd) ! Qlm_gg of density
            END IF
            nchan_0 = nsoset(lmax0)

            IF (nchan_0 > max_iso_not0) CPABORT("channels for rho0 > # max of spherical harmonics")

            NULLIFY (rrad_z, my_CG)
            my_CG => harmonics%my_CG

            !     set to zero temporary arrays
            sqrtwr = 0.0_dp
            g0_h_w = 0.0_dp
            gexp = 0.0_dp
            gsph = 0.0_dp

            sqrtwr(1:nr) = SQRT(grid_atom%wr(1:nr))
            DO l_ang = 0, lmax0
               g0_h_w(1:nr, l_ang) = g0_h(1:nr, l_ang)*grid_atom%wr(1:nr)
            END DO

            m1 = 0
            DO iset1 = 1, nset
               n1 = nsoset(lmax(iset1))
               DO ipgf1 = 1, npgf(iset1)
                  gexp(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))*sqrtwr(1:nr)
                  DO is1 = nsoset(lmin(iset1) - 1) + 1, nsoset(lmax(iset1))
                     iso = is1 + (ipgf1 - 1)*n1 + m1
                     l_ang = indso(1, is1)
                     gsph(1:nr, iso) = grid_atom%rad2l(1:nr, l_ang)*gexp(1:nr)
                  END DO ! is1
               END DO ! ipgf1
               m1 = m1 + maxso
            END DO ! iset1

            !     Distribute the atoms of this kind
            num_pe = para_env%num_pe
            mepos = para_env%mepos
            bo = get_limit(nat, num_pe, mepos)

            DO iat = bo(1), bo(2) !1,nat
               iatom = atom_list(iat)
               rho_atom => rho_atom_set(iatom)

               NULLIFY (rrad_z, vrrad_z, rrad_0, vrrad_0)
               IF (core_charge) THEN
                  rrad_z => rhoz_set(ikind)%r_coef ! for density
               END IF
               IF (my_core_2nd) THEN
                  IF (l_2nd_local_rho) THEN
                     vrrad_z => rhoz_set_2nd(ikind)%vr_coef ! for potential
                  ELSE
                     vrrad_z => rhoz_set(ikind)%vr_coef ! for potential
                  END IF
               END IF
               rrad_0 => rho0_atom_set(iatom)%rho0_rad_h%r_coef ! for density
               vrrad_0 => rho0_atom_set(iatom)%vrho0_rad_h%r_coef
               IF (l_2nd_local_rho) THEN
                  rho_atom => rho_atom_set_2nd(iatom)
                  vrrad_0 => rho0_atom_set_2nd(iatom)%vrho0_rad_h%r_coef ! for potential
               END IF
               IF (my_periodic .AND. back_ch .GT. 1.E-3_dp) THEN
                  factor = -2.0_dp*pi/3.0_dp*SQRT(fourpi)*qs_charges%background
               ELSE
                  factor = 0._dp
               END IF

               CALL Vh_1c_atom_potential(rho_atom, vrrad_0, &
                                         grid_atom, my_core_2nd, vrrad_z, Vh1_h, Vh1_s, & ! core charge for potential (2nd)
                                         nchan_0, nspins, max_iso_not0, factor)

               IF (l_2nd_local_rho) rho_atom => rho_atom_set(iatom) ! rho_atom for density

               CALL Vh_1c_atom_energy(energy_hartree_1c, ecoul_1c, rho_atom, rrad_0, &
                                      grid_atom, iatom, core_charge, rrad_z, Vh1_h, Vh1_s, & ! core charge for density
                                      nchan_0, nspins, max_iso_not0)

               IF (l_2nd_local_rho) rho_atom => rho_atom_set_2nd(iatom) ! rho_atom for potential (2nd)

               CALL Vh_1c_atom_integrals(rho_atom, &  ! results (int_local_h and int_local_s) written on rho_atom_2nd
                                         ! int_local_h and int_local_s are used in update_ks_atom
                                         ! on int_local_h mixed core / non-core
                                         aVh1b_hh, aVh1b_ss, aVh1b_00, Vh1_h, Vh1_s, max_iso_not0, &
                                         max_s_harm, llmax, cg_list, cg_n_list, &
                                         nset, npgf, lmin, lmax, nsotot, maxso, nspins, nchan_0, gsph, &
                                         g0_h_w, my_CG, Qlm_gg) ! Qlm_gg for density from local_rho_set

            END DO ! iat

            DEALLOCATE (aVh1b_hh)
            DEALLOCATE (aVh1b_ss)
            DEALLOCATE (aVh1b_00)
            DEALLOCATE (Vh1_h, Vh1_s)
            DEALLOCATE (cg_list, cg_n_list)
            DEALLOCATE (gsph)
            DEALLOCATE (gexp)
            DEALLOCATE (sqrtwr, g0_h_w)

         ELSE
!===========   NO  PAW   ===============
!  This term is taken care of using the core density as in GPW
            CYCLE
         END IF ! paw
      END DO ! ikind

      CALL para_env%sum(energy_hartree_1c)

      CALL timestop(handle)

   END SUBROUTINE Vh_1c_gg_integrals

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param rho_atom ...
!> \param vrrad_0 ...
!> \param grid_atom ...
!> \param core_charge ...
!> \param vrrad_z ...
!> \param Vh1_h ...
!> \param Vh1_s ...
!> \param nchan_0 ...
!> \param nspins ...
!> \param max_iso_not0 ...
!> \param bfactor ...
! **************************************************************************************************
   SUBROUTINE Vh_1c_atom_potential(rho_atom, vrrad_0, &
                                   grid_atom, core_charge, vrrad_z, Vh1_h, Vh1_s, &
                                   nchan_0, nspins, max_iso_not0, bfactor)

      TYPE(rho_atom_type), POINTER                       :: rho_atom
      REAL(dp), DIMENSION(:, :), POINTER                 :: vrrad_0
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      LOGICAL, INTENT(IN)                                :: core_charge
      REAL(dp), DIMENSION(:), POINTER                    :: vrrad_z
      REAL(dp), DIMENSION(:, :), POINTER                 :: Vh1_h, Vh1_s
      INTEGER, INTENT(IN)                                :: nchan_0, nspins, max_iso_not0
      REAL(dp), INTENT(IN)                               :: bfactor

      INTEGER                                            :: ir, iso, ispin, nr
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: vr_h, vr_s

      nr = grid_atom%nr

      NULLIFY (vr_h, vr_s)
      CALL get_rho_atom(rho_atom=rho_atom, vrho_rad_h=vr_h, vrho_rad_s=vr_s)

      Vh1_h = 0.0_dp
      Vh1_s = 0.0_dp

      IF (core_charge) Vh1_h(:, 1) = vrrad_z(:)

      DO iso = 1, nchan_0
         Vh1_s(:, iso) = vrrad_0(:, iso)
      END DO

      DO ispin = 1, nspins
         DO iso = 1, max_iso_not0
            Vh1_h(:, iso) = Vh1_h(:, iso) + vr_h(ispin)%r_coef(:, iso)
            Vh1_s(:, iso) = Vh1_s(:, iso) + vr_s(ispin)%r_coef(:, iso)
         END DO
      END DO

      IF (bfactor /= 0._dp) THEN
         DO ir = 1, nr
            Vh1_h(ir, 1) = Vh1_h(ir, 1) + bfactor*grid_atom%rad2(ir)*grid_atom%wr(ir)
            Vh1_s(ir, 1) = Vh1_s(ir, 1) + bfactor*grid_atom%rad2(ir)*grid_atom%wr(ir)
         END DO
      END IF

   END SUBROUTINE Vh_1c_atom_potential

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param energy_hartree_1c ...
!> \param ecoul_1c ...
!> \param rho_atom ...
!> \param rrad_0 ...
!> \param grid_atom ...
!> \param iatom ...
!> \param core_charge ...
!> \param rrad_z ...
!> \param Vh1_h ...
!> \param Vh1_s ...
!> \param nchan_0 ...
!> \param nspins ...
!> \param max_iso_not0 ...
! **************************************************************************************************
   SUBROUTINE Vh_1c_atom_energy(energy_hartree_1c, ecoul_1c, rho_atom, rrad_0, &
                                grid_atom, iatom, core_charge, rrad_z, Vh1_h, Vh1_s, &
                                nchan_0, nspins, max_iso_not0)

      REAL(dp), INTENT(INOUT)                            :: energy_hartree_1c
      TYPE(ecoul_1center_type), DIMENSION(:), POINTER    :: ecoul_1c
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      REAL(dp), DIMENSION(:, :), POINTER                 :: rrad_0
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      INTEGER, INTENT(IN)                                :: iatom
      LOGICAL, INTENT(IN)                                :: core_charge
      REAL(dp), DIMENSION(:), POINTER                    :: rrad_z
      REAL(dp), DIMENSION(:, :), POINTER                 :: Vh1_h, Vh1_s
      INTEGER, INTENT(IN)                                :: nchan_0, nspins, max_iso_not0

      INTEGER                                            :: iso, ispin, nr
      REAL(dp)                                           :: ecoul_1_0, ecoul_1_h, ecoul_1_s, &
                                                            ecoul_1_z
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: r_h, r_s

      nr = grid_atom%nr

      NULLIFY (r_h, r_s)
      CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s)

      !       Calculate the contributions to Ecoul coming from Vh1_h*rhoz
      ecoul_1_z = 0.0_dp
      IF (core_charge) THEN
         ecoul_1_z = 0.5_dp*SUM(Vh1_h(:, 1)*rrad_z(:)*grid_atom%wr(:))
      END IF

      !       Calculate the contributions to Ecoul coming from  Vh1_s*rho0
      ecoul_1_0 = 0.0_dp
      DO iso = 1, nchan_0
         ecoul_1_0 = ecoul_1_0 + 0.5_dp*SUM(Vh1_s(:, iso)*rrad_0(:, iso)*grid_atom%wr(:))
      END DO

      !       Calculate the contributions to Ecoul coming from Vh1_h*rho1_h and Vh1_s*rho1_s
      ecoul_1_s = 0.0_dp
      ecoul_1_h = 0.0_dp
      DO ispin = 1, nspins
         DO iso = 1, max_iso_not0
            ecoul_1_s = ecoul_1_s + 0.5_dp*SUM(Vh1_s(:, iso)*r_s(ispin)%r_coef(:, iso)*grid_atom%wr(:))
            ecoul_1_h = ecoul_1_h + 0.5_dp*SUM(Vh1_h(:, iso)*r_h(ispin)%r_coef(:, iso)*grid_atom%wr(:))
         END DO
      END DO

      CALL set_ecoul_1c(ecoul_1c, iatom, ecoul_1_z=ecoul_1_z, ecoul_1_0=ecoul_1_0)
      CALL set_ecoul_1c(ecoul_1c=ecoul_1c, iatom=iatom, ecoul_1_h=ecoul_1_h, ecoul_1_s=ecoul_1_s)

      energy_hartree_1c = energy_hartree_1c + ecoul_1_z - ecoul_1_0
      energy_hartree_1c = energy_hartree_1c + ecoul_1_h - ecoul_1_s

   END SUBROUTINE Vh_1c_atom_energy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! **************************************************************************************************
!> \brief ...
!> \param rho_atom ...
!> \param aVh1b_hh ...
!> \param aVh1b_ss ...
!> \param aVh1b_00 ...
!> \param Vh1_h ...
!> \param Vh1_s ...
!> \param max_iso_not0 ...
!> \param max_s_harm ...
!> \param llmax ...
!> \param cg_list ...
!> \param cg_n_list ...
!> \param nset ...
!> \param npgf ...
!> \param lmin ...
!> \param lmax ...
!> \param nsotot ...
!> \param maxso ...
!> \param nspins ...
!> \param nchan_0 ...
!> \param gsph ...
!> \param g0_h_w ...
!> \param my_CG ...
!> \param Qlm_gg ...
! **************************************************************************************************
   SUBROUTINE Vh_1c_atom_integrals(rho_atom, &
                                   aVh1b_hh, aVh1b_ss, aVh1b_00, Vh1_h, Vh1_s, max_iso_not0, &
                                   max_s_harm, llmax, cg_list, cg_n_list, &
                                   nset, npgf, lmin, lmax, nsotot, maxso, nspins, nchan_0, gsph, &
                                   g0_h_w, my_CG, Qlm_gg)

      TYPE(rho_atom_type), POINTER                       :: rho_atom
      REAL(dp), DIMENSION(:, :)                          :: aVh1b_hh, aVh1b_ss, aVh1b_00
      REAL(dp), DIMENSION(:, :), POINTER                 :: Vh1_h, Vh1_s
      INTEGER, INTENT(IN)                                :: max_iso_not0, max_s_harm, llmax
      INTEGER, DIMENSION(:, :, :)                        :: cg_list
      INTEGER, DIMENSION(:)                              :: cg_n_list
      INTEGER, INTENT(IN)                                :: nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, lmin, lmax
      INTEGER, INTENT(IN)                                :: nsotot, maxso, nspins, nchan_0
      REAL(dp), DIMENSION(:, :), POINTER                 :: gsph
      REAL(dp), DIMENSION(:, 0:)                         :: g0_h_w
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG, Qlm_gg

      INTEGER                                            :: icg, ipgf1, ipgf2, ir, is1, is2, iset1, &
                                                            iset2, iso, iso1, iso2, ispin, l_ang, &
                                                            m1, m2, max_iso_not0_local, n1, n2, nr
      REAL(dp)                                           :: gVg_0, gVg_h, gVg_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_local_h, int_local_s

      NULLIFY (int_local_h, int_local_s)
      CALL get_rho_atom(rho_atom=rho_atom, &
                        ga_Vlocal_gb_h=int_local_h, &
                        ga_Vlocal_gb_s=int_local_s)

      !       Calculate the integrals of the potential with 2 primitives
      aVh1b_hh = 0.0_dp
      aVh1b_ss = 0.0_dp
      aVh1b_00 = 0.0_dp

      nr = SIZE(gsph, 1)

      m1 = 0
      DO iset1 = 1, nset
         m2 = 0
         DO iset2 = 1, nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, llmax, cg_list, cg_n_list, max_iso_not0_local)

            n1 = nsoset(lmax(iset1))
            DO ipgf1 = 1, npgf(iset1)
               n2 = nsoset(lmax(iset2))
               DO ipgf2 = 1, npgf(iset2)
                  !               with contributions to  V1_s*rho0
                  DO iso = 1, nchan_0
                     l_ang = indso(1, iso)
                     gVg_0 = SUM(Vh1_s(:, iso)*g0_h_w(:, l_ang))
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2
                        gVg_h = 0.0_dp
                        gVg_s = 0.0_dp

                        DO ir = 1, nr
                           gVg_h = gVg_h + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_h(ir, iso)
                           gVg_s = gVg_s + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_s(ir, iso)
                        END DO ! ir

                        aVh1b_hh(iso1, iso2) = aVh1b_hh(iso1, iso2) + gVg_h*my_CG(is1, is2, iso)
                        aVh1b_ss(iso1, iso2) = aVh1b_ss(iso1, iso2) + gVg_s*my_CG(is1, is2, iso)
                        aVh1b_00(iso1, iso2) = aVh1b_00(iso1, iso2) + gVg_0*Qlm_gg(iso1, iso2, iso)

                     END DO !icg
                  END DO ! iso
                  !               without contributions to  V1_s*rho0
                  DO iso = nchan_0 + 1, max_iso_not0
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2
                        gVg_h = 0.0_dp
                        gVg_s = 0.0_dp

                        DO ir = 1, nr
                           gVg_h = gVg_h + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_h(ir, iso)
                           gVg_s = gVg_s + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_s(ir, iso)
                        END DO ! ir

                        aVh1b_hh(iso1, iso2) = aVh1b_hh(iso1, iso2) + gVg_h*my_CG(is1, is2, iso)
                        aVh1b_ss(iso1, iso2) = aVh1b_ss(iso1, iso2) + gVg_s*my_CG(is1, is2, iso)

                     END DO !icg
                  END DO ! iso
               END DO ! ipgf2
            END DO ! ipgf1
            m2 = m2 + maxso
         END DO ! iset2
         m1 = m1 + maxso
      END DO !iset1
      DO ispin = 1, nspins
         CALL daxpy(nsotot*nsotot, 1.0_dp, aVh1b_hh, 1, int_local_h(ispin)%r_coef, 1)
         CALL daxpy(nsotot*nsotot, 1.0_dp, aVh1b_ss, 1, int_local_s(ispin)%r_coef, 1)
         CALL daxpy(nsotot*nsotot, -1.0_dp, aVh1b_00, 1, int_local_h(ispin)%r_coef, 1)
         CALL daxpy(nsotot*nsotot, -1.0_dp, aVh1b_00, 1, int_local_s(ispin)%r_coef, 1)
      END DO ! ispin

   END SUBROUTINE Vh_1c_atom_integrals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE hartree_local_methods

