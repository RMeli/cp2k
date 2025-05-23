!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_tb
   USE bibliography,                    ONLY: Elstner1998,&
                                              Grimme2017,&
                                              Hu2007,&
                                              Porezag1995,&
                                              Seifert1996,&
                                              Zhechkov2005
   USE eeq_input,                       ONLY: create_eeq_control_section
   USE input_constants,                 ONLY: dispersion_d2,&
                                              dispersion_d3,&
                                              dispersion_d3bj,&
                                              dispersion_uff,&
                                              gfn1xtb,&
                                              gfn2xtb,&
                                              ipea1xtb,&
                                              slater
   USE input_cp2k_mm,                   ONLY: create_GENPOT_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tb'

   PUBLIC :: create_dftb_control_section, create_xtb_control_section

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_dftb_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DFTB", &
                          description="Parameters needed to set up the DFTB methods", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=(/Porezag1995, Seifert1996, Elstner1998, Zhechkov2005/))

      NULLIFY (subsection)
      CALL create_dftb_parameter_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="self_consistent", &
                          description="Use self-consistent method", &
                          citations=(/Elstner1998/), &
                          usage="SELF_CONSISTENT", default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="orthogonal_basis", &
                          description="Assume orthogonal basis set", &
                          usage="ORTHOGONAL_BASIS", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="do_ewald", &
                          description="Use Ewald type method instead of direct sum for Coulomb interaction", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="dispersion", &
                          description="Use dispersion correction", &
                          citations=(/Zhechkov2005/), lone_keyword_l_val=.TRUE., &
                          usage="DISPERSION", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIAGONAL_DFTB3", &
                          description="Use a diagonal version of the 3rd order energy correction (DFTB3) ", &
                          lone_keyword_l_val=.TRUE., &
                          usage="DIAGONAL_DFTB3", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HB_SR_GAMMA", &
                          description="Uses a modified version for the GAMMA within the SCC-DFTB scheme, "// &
                          "specifically tuned for hydrogen bonds.", &
                          citations=(/Hu2007/), lone_keyword_l_val=.TRUE., &
                          usage="HB_SR_GAMMA", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_disp", &
                          description="Define accuracy of dispersion interaction", &
                          usage="EPS_DISP", default_r_val=0.0001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dftb_control_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xTB", &
                          description="Parameters needed to set up the xTB methods", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=(/GRIMME2017/))

      NULLIFY (subsection)
      CALL create_xtb_parameter_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xtb_nonbonded_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eeq_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xtb_tblite_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GFN_TYPE", &
                          description="Which GFN xTB method should be used.", &
                          usage="GFN_TYPE 1", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_EWALD", &
                          description="Use Ewald type method instead of direct sum for Coulomb interaction", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STO_NG", &
                          description="Provides the order of the Slater orbital expansion in GTOs.", &
                          usage="STO_NG 3", default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HYDROGEN_STO_NG", &
                          description="Number of GTOs for Hydrogen basis expansion.", &
                          usage="HYDROGEN_STO_NG 3", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_HALOGEN_CORRECTION", &
                          description="Use XB interaction term", &
                          usage="USE_HALOGEN_CORRECTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_NONBONDED", &
                          description="Controls the computation of real-space "// &
                          "(short-range) nonbonded interactions as correction to xTB.", &
                          usage="DO_NONBONDED T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VDW_POTENTIAL", &
                          description="vdW potential to be used: NONE, DFTD3, DFTD4. "// &
                          "Defaults: DFTD3(gfn1), DFTD4(gfn0, gfn2).", &
                          usage="VDW_POTENTIAL type", default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_INTERACTION", &
                          description="Use Coulomb interaction terms (electrostatics + TB3); for debug only", &
                          usage="COULOMB_INTERACTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_LR", &
                          description="Use Coulomb LR (1/r) interaction terms; for debug only", &
                          usage="COULOMB_LR T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TB3_INTERACTION", &
                          description="Use TB3 interaction terms; for debug only", &
                          usage="TB3_INTERACTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_ATOMIC_CHARGES", &
                          description="Stop calculation if atomic charges are outside chemical range.", &
                          usage="CHECK_ATOMIC_CHARGES T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VARIATIONAL_DIPOLE", &
                          description="gfn0-xTB use dipole definition from energy derivative.", &
                          usage="VARIATIONAL_DIPOLE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PAIRPOTENTIAL", &
                          description="Accuracy for the repulsive pair potential.", &
                          usage="EPS_PAIRPOTENTIAL 1.0E-8", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EN_SHIFT_TYPE", &
                          description="Shift function for electronegativity in EEQ method. "// &
                          "[Select/Molecule/Crystal] Default Select from periodicity.", &
                          usage="EN_SHIFT_TYPE [Select/Molecule/Crystal]", &
                          n_var=1, type_of_var=char_t, default_c_val="Molecule")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_control_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_dftb_parameter_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="PARAMETER", &
                          description="Information on where to find DFTB parameters", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SK_FILE", &
                          description="Define parameter file for atom pair", &
                          usage="SK_FILE a1 a2 filename", &
                          n_var=3, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_PATH", &
                          description="Specify the directory with the DFTB parameter files. "// &
                          "Used in combination with the filenames specified in the file "// &
                          "given in PARAM_FILE_NAME.", usage="PARAM_FILE_PATH pathname", &
                          n_var=1, type_of_var=char_t, default_c_val="./")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_NAME", &
                          description="Specify file that contains the names of "// &
                          "Slater-Koster tables: A plain text file, each line has the "// &
                          'format "ATOM1 ATOM2 filename.spl".', &
                          usage="PARAM_FILE_NAME filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_TYPE", &
                          description="Use dispersion correction of the specified type."// &
                          " Dispersion correction has to be switched on in the DFTB section.", &
                          usage="DISPERSION_TYPE (UFF|D3|D3(BJ)|D2)", &
                          enum_c_vals=s2a("UFF", "D3", "D3(BJ)", "D2"), &
                          enum_i_vals=(/dispersion_uff, dispersion_d3, dispersion_d3bj, dispersion_d2/), &
                          enum_desc=s2a("Uses the UFF force field for a pair potential dispersion correction.", &
                                        "Uses the Grimme D3 method (simplified) for a pair potential dispersion correction.", &
                                        "Uses the Grimme D3 method (simplified) with Becke-Johnson attenuation.", &
                                        "Uses the Grimme D2 method for pair potential dispersion correction."), &
                          default_i_val=dispersion_uff)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UFF_FORCE_FIELD", &
                          description="Name of file with UFF parameters that will be used "// &
                          "for the dispersion correction. Needs to be specified when "// &
                          "DISPERSION==.TRUE., otherwise cp2k crashes with a Segmentation "// &
                          "Fault.", usage="UFF_FORCE_FIELD filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion "// &
                          "parameters for the D3 method", &
                          usage="DISPERSION_PARAMETER_FILE filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_RADIUS", &
                          description="Define radius of dispersion interaction", &
                          usage="DISPERSION_RADIUS", default_r_val=15._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COORDINATION_CUTOFF", &
                          description="Define cutoff for coordination number calculation", &
                          usage="COORDINATION_CUTOFF", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3_SCALING", &
                          description="Scaling parameters (s6,sr6,s8) for the D3 dispersion method,", &
                          usage="D3_SCALING 1.0 1.0 1.0", n_var=3, default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_SCALING", &
                          description="Scaling parameters (s6,a1,s8,a2) for the D3(BJ) dispersion method,", &
                          usage="D3BJ_SCALING 1.0 1.0 1.0 1.0", n_var=4, &
                          default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D2_SCALING", &
                          description="Scaling parameter for the D2 dispersion method,", &
                          usage="D2_SCALING 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D2_EXP_PRE", &
                          description="Exp prefactor for damping for the D2 dispersion method,", &
                          usage="D2_EXP_PRE 2.0", default_r_val=2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HB_SR_PARAM", &
                          description="Uses a modified version for the GAMMA within the SCC-DFTB scheme, "// &
                          "specifically tuned for hydrogen bonds. Specify the exponent used in the exponential.", &
                          usage="HB_SR_PARAM {real}", default_r_val=4.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dftb_parameter_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_parameter_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="PARAMETER", &
                          description="Information on and where to find xTB parameters", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_PATH", &
                          description="Specify the directory with the xTB parameter file. ", &
                          usage="PARAM_FILE_PATH pathname", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_NAME", &
                          description="Specify file that contains all xTB default parameters. ", &
                          usage="PARAM_FILE_NAME filename", &
                          n_var=1, type_of_var=char_t, default_c_val="xTB_parameters")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion "// &
                          "parameters for the D3 method", &
                          usage="DISPERSION_PARAMETER_FILE filename", &
                          n_var=1, type_of_var=char_t, default_c_val="dftd3.dat")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_RADIUS", &
                          description="Define radius of dispersion interaction", &
                          usage="DISPERSION_RADIUS", default_r_val=15._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COORDINATION_CUTOFF", &
                          description="Define cutoff for coordination number calculation", &
                          usage="COORDINATION_CUTOFF", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_SCALING", &
                          description="Scaling parameters (s6,s8) for the D3 dispersion method.", &
                          usage="D3BJ_SCALING 1.0 2.4", n_var=2, default_r_vals=(/1.0_dp, 2.4_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_PARAM", &
                          description="Becke-Johnson parameters (a1, a2 for the D3 dispersion method.", &
                          usage="D3BJ_PARAM 0.63 5.0", n_var=2, default_r_vals=(/0.63_dp, 5.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HUCKEL_CONSTANTS", &
                          description="Huckel parameters (s, p, d, sp, 2sH).", &
                          usage="HUCKEL_CONSTANTS 1.85 2.25 2.00 2.08 2.85", n_var=5, &
                          default_r_vals=(/1.85_dp, 2.25_dp, 2.00_dp, 2.08_dp, 2.85_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_CONSTANTS", &
                          description="Scaling parameters for Coulomb interactions (electrons, nuclei).", &
                          usage="COULOMB_CONSTANTS 2.00 1.50", n_var=2, &
                          default_r_vals=(/2.00_dp, 1.50_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CN_CONSTANTS", &
                          description="Scaling parameters for Coordination number correction term.", &
                          usage="CN_CONSTANTS 0.006 -0.003 -0.005", n_var=3, &
                          default_r_vals=(/0.006_dp, -0.003_dp, -0.005_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EN_CONSTANTS", &
                          description="Scaling parameters for electronegativity correction term.", &
                          usage="EN_CONSTANTS -0.007 0.000 0.000", n_var=3, &
                          default_r_vals=(/-0.007_dp, 0.000_dp, 0.000_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BEN_CONSTANT", &
                          description="Scaling parameter for electronegativity correction term.", &
                          usage="BEN_CONSTANT 4.0", n_var=1, &
                          default_r_val=4.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENSCALE", &
                          description="Scaling parameter repulsive energy (dEN in exponential).", &
                          usage="ENSCALE 0.01", n_var=1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HALOGEN_BINDING", &
                          description="Scaling parameters for electronegativity correction term.", &
                          usage="HALOGEN_BINDING 1.30 0.44", n_var=2, default_r_vals=(/1.30_dp, 0.44_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KAB_PARAM", &
                          description="Specifies the specific Kab value for types A and B.", &
                          usage="KAB_PARAM kind1 kind2 value ", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XB_RADIUS", &
                          description="Specifies the radius [Bohr] of the XB pair interaction in xTB.", &
                          usage="XB_RADIUS 20.0 ", repeats=.FALSE., &
                          n_var=1, default_r_val=20.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_SR_CUT", &
                          description="Maximum range of short range part of Coulomb interaction.", &
                          usage="COULOMB_SR_CUT 20.0 ", repeats=.FALSE., &
                          n_var=1, default_r_val=20.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_SR_EPS", &
                          description="Cutoff for short range part of Coulomb interaction.", &
                          usage="COULOMB_SR_EPS 1.E-3 ", repeats=.FALSE., &
                          n_var=1, default_r_val=1.0E-03_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SRB_PARAMETER", &
                          description="SRB parameters (ksrb, esrb, gscal, c1, c2, shift).", &
                          usage="SRB_PARAMETER -0.0129 3.48 0.51 -1.71 2.11 0.0537", n_var=6, &
                          default_r_vals=(/-0.0129_dp, 3.4847_dp, 0.5097_dp, &
                                           -1.70549806_dp, 2.10878369_dp, 0.0537_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_parameter_section
! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_nonbonded_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="NONBONDED", &
                          description="This section specifies the input parameters for NON-BONDED interactions.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (subsection)

      CALL create_GENPOT_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_nonbonded_section
! **************************************************************************************************
!> \brief Creates the &TBLITE section
!> \param section the section to create
!> \author JVP
! **************************************************************************************************
   SUBROUTINE create_xtb_tblite_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TBLITE", &
                          description="Section used to specify options for an xTB computation using tblite.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="activates the execution via tblite", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Selection of the method used in tblite.", &
                          usage="METHOD (GFN1|GFN2|IPEA1)", &
                          enum_c_vals=s2a("GFN1", "GFN2", "IPEA1"), &
                          enum_i_vals=(/gfn1xtb, gfn2xtb, ipea1xtb/), &
                          enum_desc=s2a("Uses the GFN1-XTB method by Grimme.", &
                                        "Uses the GFN2-XTB method by Grimme.", &
                                        "Uses the IEPEA1 method by Grimme."), &
                          default_i_val=gfn2xtb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_tblite_section

END MODULE input_cp2k_tb
