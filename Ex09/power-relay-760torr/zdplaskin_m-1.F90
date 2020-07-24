!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Thu May  7 19:22:28 2020
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E C O 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 42, species_electrons = 1, species_length = 7, reactions_max = 503, reactions_length = 49
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 23, bolsig_species_length = 7, bolsig_rates_max = 147 
  character(*), parameter, private          :: bolsigfile = "CO2SPLIT2C_BOLSIGDBEE.DAT"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  /-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1, 1, 1, 1,-1,-1,-1,-1/
  data species_name(1:species_max) &
  /"E      ","CO2    ","CO2(VA)","CO2(VB)","CO2(V1)","CO2(VD)","CO2(E1)","CO2(E2)","CO     ","CO(V1) ","CO(E1) ","CO(E2) ",&
   "CO(E3) ","CO(E4) ","O2     ","O2(V1) ","O2(V2) ","O2(V3) ","O2(E1) ","O2(E2) ","C2O    ","C      ","C2     ","O3     ",&
   "O      ","CO2^+  ","CO4^+  ","CO^+   ","C2O2^+ ","C2O3^+ ","C2O4^+ ","C^+    ","C2^+   ","CO3^-  ","CO4^-  ","O2^+   ",&
   "O^+    ","O4^+   ","O^-    ","O2^-   ","O3^-   ","O4^-   "/
  data reaction_sign(1:36) &
  /"bolsig:CO2->CO2^+                                ","bolsig:CO2->CO^+                                 ",&
   "bolsig:CO2->C^+                                  ","bolsig:CO2->O^+                                  ",&
   "E+CO2=>O2^++C+E+E                                ","bolsig:CO2->CO+O^-                               ",&
   "bolsig:CO2->CO+O                                 ","bolsig:CO2->CO2(VA)                              ",&
   "bolsig:CO2->CO2(VB)                              ","bolsig:CO2->CO2(V1)                              ",&
   "bolsig:CO2->CO2(VD)                              ","bolsig:CO2->CO2(E1)                              ",&
   "bolsig:CO2->CO2(E2)                              ","bolsig:CO2(E1)->CO2                              ",&
   "bolsig:CO2(E2)->CO2                              ","bolsig:CO2(VA)->CO+O                             ",&
   "bolsig:CO2(VB)->CO+O                             ","bolsig:CO2(VD)->CO+O                             ",&
   "bolsig:CO2(V1)->CO+O                             ","bolsig:CO2(VA)->CO2^+                            ",&
   "bolsig:CO2(VB)->CO2^+                            ","bolsig:CO2(VD)->CO2^+                            ",&
   "bolsig:CO2(V1)->CO2^+                            ","bolsig:CO2(VA)->CO^+                             ",&
   "bolsig:CO2(VB)->CO^+                             ","bolsig:CO2(VD)->CO^+                             ",&
   "bolsig:CO2(V1)->CO^+                             ","bolsig:CO2(VA)->C^+                              ",&
   "bolsig:CO2(VB)->C^+                              ","bolsig:CO2(VD)->C^+                              ",&
   "bolsig:CO2(V1)->C^+                              ","bolsig:CO2(VA)->O^+                              ",&
   "bolsig:CO2(VB)->O^+                              ","bolsig:CO2(VD)->O^+                              ",&
   "bolsig:CO2(V1)->O^+                              ","bolsig:CO2(VA)->CO+O^-                           "/
  data reaction_sign(37:72) &
  /"bolsig:CO2(VB)->CO+O^-                           ","bolsig:CO2(VD)->CO+O^-                           ",&
   "bolsig:CO2(V1)->CO+O^-                           ","bolsig:CO2(VA)->CO2(E1)                          ",&
   "bolsig:CO2(VB)->CO2(E1)                          ","bolsig:CO2(VD)->CO2(E1)                          ",&
   "bolsig:CO2(V1)->CO2(E1)                          ","bolsig:CO2(VA)->CO2(E2)                          ",&
   "bolsig:CO2(VB)->CO2(E2)                          ","bolsig:CO2(VD)->CO2(E2)                          ",&
   "bolsig:CO2(V1)->CO2(E2)                          ","E+CO2(VA)=>O2^++C+E+E                            ",&
   "E+CO2(VB)=>O2^++C+E+E                            ","E+CO2(VD)=>O2^++C+E+E                            ",&
   "E+CO2(V1)=>O2^++C+E+E                            ","bolsig:CO2(VA)->CO2                              ",&
   "bolsig:CO2(VB)->CO2                              ","bolsig:CO2(VD)->CO2                              ",&
   "bolsig:CO2(V1)->CO2                              ","bolsig:CO2(E1)->CO+O                             ",&
   "bolsig:CO2(E2)->CO+O                             ","bolsig:CO2(E1)->CO2^+                            ",&
   "bolsig:CO2(E2)->CO2^+                            ","bolsig:CO2(E1)->CO^+                             ",&
   "bolsig:CO2(E2)->CO^+                             ","bolsig:CO2(E1)->C^+                              ",&
   "bolsig:CO2(E2)->C^+                              ","bolsig:CO2(E1)->O^+                              ",&
   "bolsig:CO2(E2)->O^+                              ","E+CO2(E1)=>O2^++C+E+E                            ",&
   "E+CO2(E2)=>O2^++C+E+E                            ","bolsig:CO->CO^+                                  ",&
   "bolsig:CO->C^+                                   ","bolsig:CO->O^+                                   ",&
   "bolsig:CO->O^-                                   ","bolsig:CO->C+O                                   "/
  data reaction_sign(73:108) &
  /"bolsig:CO->CO(V1)                                ","bolsig:CO->CO(E1)                                ",&
   "bolsig:CO->CO(E2)                                ","bolsig:CO->CO(E3)                                ",&
   "bolsig:CO->CO(E4)                                ","bolsig:CO(E1)->CO                                ",&
   "bolsig:CO(E2)->CO                                ","bolsig:CO(E3)->CO                                ",&
   "bolsig:CO(E4)->CO                                ","bolsig:CO(V1)->C+O                               ",&
   "bolsig:CO(V1)->CO^+                              ","bolsig:CO(V1)->C^+                               ",&
   "bolsig:CO(V1)->O^+                               ","bolsig:CO(V1)->O^-                               ",&
   "bolsig:CO(V1)->CO(E1)                            ","bolsig:CO(V1)->CO(E2)                            ",&
   "bolsig:CO(V1)->CO(E3)                            ","bolsig:CO(V1)->CO(E4)                            ",&
   "bolsig:CO(V1)->CO                                ","bolsig:CO(E1)->C+O                               ",&
   "bolsig:CO(E2)->C+O                               ","bolsig:CO(E3)->C+O                               ",&
   "bolsig:CO(E4)->C+O                               ","bolsig:CO(E1)->CO^+                              ",&
   "bolsig:CO(E2)->CO^+                              ","bolsig:CO(E3)->CO^+                              ",&
   "bolsig:CO(E4)->CO^+                              ","bolsig:CO(E1)->C^+                               ",&
   "bolsig:CO(E2)->C^+                               ","bolsig:CO(E3)->C^+                               ",&
   "bolsig:CO(E4)->C^+                               ","bolsig:CO(E1)->O^+                               ",&
   "bolsig:CO(E2)->O^+                               ","bolsig:CO(E3)->O^+                               ",&
   "bolsig:CO(E4)->O^+                               ","bolsig:C->C^+                                    "/
  data reaction_sign(109:144) &
  /"bolsig:C2->C+C                                   ","bolsig:C2->C2^+                                  ",&
   "bolsig:O2->O+O                                   ","bolsig:O2->O2^+                                  ",&
   "bolsig:O2->O^+                                   ","bolsig:O2->O^-                                   ",&
   "bolsig:O2->O2(V1)                                ","bolsig:O2->O2(V2)                                ",&
   "bolsig:O2->O2(V3)                                ","bolsig:O2->O2(E1)                                ",&
   "bolsig:O2->O2(E2)                                ","bolsig:O2(E1)->O2                                ",&
   "bolsig:O2(E2)->O2                                ","O2(E1)=>O2                                       ",&
   "O2(E2)=>O2                                       ","bolsig:O2->O2(E1)                                ",&
   "bolsig:O2->O2(E1)                                ","bolsig:O2->O2(E1)                                ",&
   "bolsig:O2->O2(E2)                                ","bolsig:O2->O2(E2)                                ",&
   "bolsig:O2->O2(E2)                                ","bolsig:O2->O2^+                                  ",&
   "bolsig:O2->O2^+                                  ","bolsig:O2->O2^+                                  ",&
   "bolsig:O2->O^+                                   ","bolsig:O2->O^+                                   ",&
   "bolsig:O2->O^+                                   ","bolsig:O2->O+O                                   ",&
   "bolsig:O2->O+O                                   ","bolsig:O2->O+O                                   ",&
   "bolsig:O2->O^-                                   ","bolsig:O2->O^-                                   ",&
   "bolsig:O2->O^-                                   ","bolsig:O2(V1)->O2                                ",&
   "bolsig:O2(V2)->O2                                ","bolsig:O2(V3)->O2                                "/
  data reaction_sign(145:180) &
  /"bolsig:O2(E1)->O+O                               ","bolsig:O2(E2)->O+O                               ",&
   "bolsig:O2(E1)->O2^+                              ","bolsig:O2(E2)->O2^+                              ",&
   "bolsig:O2(E1)->O^+                               ","bolsig:O2(E2)->O^+                               ",&
   "bolsig:O3->O2+O                                  ","bolsig:O3->O2^+                                  ",&
   "bolsig:O3->O^+                                   ","bolsig:O3->O^-                                   ",&
   "bolsig:O3->O2^-                                  ","bolsig:O->O^+                                    ",&
   "E+CO2^+=>CO+O                                    ","E+CO2^+=>C+O2                                    ",&
   "E+CO4^+=>CO2+O2                                  ","E+CO^+=>C+O                                      ",&
   "E+C2O2^+=>CO+CO                                  ","E+C2O3^+=>CO2+CO                                 ",&
   "E+C2O4^+=>CO2+CO2                                ","E+C2^+=>C+C                                      ",&
   "E+O2+ANY_NEUTRAL=>O2^-+ANY_NEUTRAL               ","E+O3+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL               ",&
   "E+O+ANY_NEUTRAL=>O^-+ANY_NEUTRAL                 ","E+O2^++ANY_NEUTRAL=>O2+ANY_NEUTRAL               ",&
   "E+O2^+=>O+O                                      ","E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL                 ",&
   "E+O4^+=>O2+O2                                    ","O+CO2=>CO+O2                                     ",&
   "C+CO2=>CO+CO                                     ","CO2(VA)+O=>CO+O2                                 ",&
   "CO2(VB)+O=>CO+O2                                 ","CO2(V1)+O=>CO+O2                                 ",&
   "CO2(VD)+O=>CO+O2                                 ","CO2(VA)+C=>2CO                                   ",&
   "CO2(VB)+C=>2CO                                   ","CO2(V1)+C=>2CO                                   "/
  data reaction_sign(181:216) &
  /"CO2(VD)+C=>2CO                                   ","O+CO+ANY_NEUTRAL=>CO2+ANY_NEUTRAL                ",&
   "O3+CO=>CO2+O2                                    ","C+CO+ANY_NEUTRAL=>C2O+ANY_NEUTRAL                ",&
   "O2+CO=>CO2+O                                     ","CO(V1)+O+ANY_NEUTRAL=>CO2+ANY_NEUTRAL            ",&
   "O2(V1)+CO=>CO2+O                                 ","O2(V2)+CO=>CO2+O                                 ",&
   "O2(V3)+CO=>CO2+O                                 ","O2+CO(V1)=>CO2+O                                 ",&
   "O2(V1)+CO(V1)=>CO2+O                             ","O2(V2)+CO(V1)=>CO2+O                             ",&
   "O2(V3)+CO(V1)=>CO2+O                             ","CO(V1)+O3=>CO2+O2                                ",&
   "CO(V1)+C+ANY_NEUTRAL=>C2O+ANY_NEUTRAL            ","O2+C=>CO+O                                       ",&
   "O2+C2O=>CO2+CO                                   ","O+O2+ANY_NEUTRAL=>O3+ANY_NEUTRAL                 ",&
   "O2(V1)+C=>CO+O                                   ","O2(V2)+C=>CO+O                                   ",&
   "O2(V3)+C=>CO+O                                   ","O2(V1)+O+ANY_NEUTRAL=>O3+ANY_NEUTRAL             ",&
   "O2(V2)+O+ANY_NEUTRAL=>O3+ANY_NEUTRAL             ","O2(V3)+O+ANY_NEUTRAL=>O3+ANY_NEUTRAL             ",&
   "O+C+ANY_NEUTRAL=>CO+ANY_NEUTRAL                  ","O+C2O=>CO+CO                                     ",&
   "O+O3=>O2+O2                                      ","O3+ANY_NEUTRAL=>O2+O+ANY_NEUTRAL                 ",&
   "O+O+ANY_NEUTRAL=>O2+ANY_NEUTRAL                  ","O2^++CO2+ANY_NEUTRAL=>CO4^++ANY_NEUTRAL          ",&
   "O^++CO2=>O2^++CO                                 ","O^++CO2=>CO2^++O                                 ",&
   "C^++CO2=>CO^++CO                                 ","CO^++CO2=>CO2^++CO                               ",&
   "O^-+CO2+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL           ","O2^-+CO2+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL          "/
  data reaction_sign(217:252) &
  /"O3^-+CO2=>O2+CO3^-                               ","O4^-+CO2=>CO4^-+O2                               ",&
   "CO2^++CO2+ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL        ","CO2(VA)+O^+=>O2^++CO                             ",&
   "CO2(VB)+O^+=>O2^++CO                             ","CO2(V1)+O^+=>O2^++CO                             ",&
   "CO2(VD)+O^+=>O2^++CO                             ","CO2(VA)+O^+=>CO2^++O                             ",&
   "CO2(VB)+O^+=>CO2^++O                             ","CO2(V1)+O^+=>CO2^++O                             ",&
   "CO2(VD)+O^+=>CO2^++O                             ","CO2(VA)+C^+=>CO^++CO                             ",&
   "CO2(VB)+C^+=>CO^++CO                             ","CO2(V1)+C^+=>CO^++CO                             ",&
   "CO2(VD)+C^+=>CO^++CO                             ","CO2(VA)+CO^+=>CO2^++CO                           ",&
   "CO2(VB)+CO^+=>CO2^++CO                           ","CO2(V1)+CO^+=>CO2^++CO                           ",&
   "CO2(VD)+CO^+=>CO2^++CO                           ","CO2(VA)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ",&
   "CO2(VB)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ","CO2(V1)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ",&
   "CO2(VD)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ","CO2(VA)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ",&
   "CO2(VB)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ","CO2(V1)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ",&
   "CO2(VD)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ","CO2(VA)+O3^-=>O2+CO3^-                           ",&
   "CO2(VB)+O3^-=>O2+CO3^-                           ","CO2(V1)+O3^-=>O2+CO3^-                           ",&
   "CO2(VD)+O3^-=>O2+CO3^-                           ","CO2(VA)+O4^-=>CO4^-+O2                           ",&
   "CO2(VB)+O4^-=>CO4^-+O2                           ","CO2(V1)+O4^-=>CO4^-+O2                           ",&
   "CO2(VD)+O4^-=>CO4^-+O2                           ","CO2(VA)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    "/
  data reaction_sign(253:288) &
  /"CO2(VB)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    ","CO2(V1)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    ",&
   "CO2(VD)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    ","CO2(VA)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ",&
   "CO2(VB)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ","CO2(V1)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ",&
   "CO2(VD)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ","CO2(E1)+O^+=>O2^++CO                             ",&
   "CO2(E2)+O^+=>O2^++CO                             ","CO2(E1)+O^+=>CO2^++O                             ",&
   "CO2(E2)+O^+=>CO2^++O                             ","CO2(E1)+C^+=>CO^++CO                             ",&
   "CO2(E2)+C^+=>CO^++CO                             ","CO2(E1)+CO^+=>CO2^++CO                           ",&
   "CO2(E2)+CO^+=>CO2^++CO                           ","CO2(E1)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ",&
   "CO2(E2)+O^-+ANY_NEUTRAL=>CO3^-+ANY_NEUTRAL       ","CO2(E1)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ",&
   "CO2(E2)+O2^-+ANY_NEUTRAL=>CO4^-+ANY_NEUTRAL      ","CO2(E1)+O3^-=>O2+CO3^-                           ",&
   "CO2(E2)+O3^-=>O2+CO3^-                           ","CO2(E1)+O4^-=>CO4^-+O2                           ",&
   "CO2(E2)+O4^-=>CO4^-+O2                           ","CO2(E1)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    ",&
   "CO2(E2)+CO2^++ANY_NEUTRAL=>C2O4^++ANY_NEUTRAL    ","CO2(E1)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ",&
   "CO2(E2)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ","CO2(E1)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ",&
   "CO2(E2)+O2^++ANY_NEUTRAL=>CO4^++ANY_NEUTRAL      ","O^++CO=>CO^++O                                   ",&
   "O^-+CO=>CO2+E                                    ","CO3^-+CO=>2CO2+E                                 ",&
   "C2O3^++CO=>CO2+C2O2^+                            ","C2O4^++CO=>C2O3^++CO2                            ",&
   "C2O3^++CO+ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL    ","C2O4^++CO+ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL    "/
  data reaction_sign(289:324) &
  /"C^++CO=>CO^++C                                   ","CO(V1)+O^+=>CO^++O                               ",&
   "CO(V1)+O^-=>CO2+E                                ","CO(V1)+CO3^-=>CO2+CO2+E                          ",&
   "CO(V1)+C2O3^+=>CO2+C2O2^+                        ","CO(V1)+C2O4^+=>C2O3^++CO2                        ",&
   "CO(V1)+C2O3^++ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL","CO(V1)+C2O4^++ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL",&
   "CO(V1)+C^+=>CO^++C                               ","CO(E1)+O^-=>CO2+E                                ",&
   "CO(E2)+O^-=>CO2+E                                ","CO(E3)+O^-=>CO2+E                                ",&
   "CO(E4)+O^-=>CO2+E                                ","CO(E1)+CO3^-=>CO2+CO2+E                          ",&
   "CO(E2)+CO3^-=>CO2+CO2+E                          ","CO(E3)+CO3^-=>CO2+CO2+E                          ",&
   "CO(E4)+CO3^-=>CO2+CO2+E                          ","CO(E1)+C2O3^+=>CO2+C2O2^+                        ",&
   "CO(E2)+C2O3^+=>CO2+C2O2^+                        ","CO(E3)+C2O3^+=>CO2+C2O2^+                        ",&
   "CO(E4)+C2O3^+=>CO2+C2O2^+                        ","CO(E1)+C2O4^+=>C2O3^++CO2                        ",&
   "CO(E2)+C2O4^+=>C2O3^++CO2                        ","CO(E3)+C2O4^+=>C2O3^++CO2                        ",&
   "CO(E4)+C2O4^+=>C2O3^++CO2                        ","CO(E1)+C2O3^++ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL",&
   "CO(E2)+C2O3^++ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL","CO(E3)+C2O3^++ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL",&
   "CO(E4)+C2O3^++ANY_NEUTRAL=>C2O2^++CO2+ANY_NEUTRAL","CO(E1)+C2O4^++ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL",&
   "CO(E2)+C2O4^++ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL","CO(E3)+C2O4^++ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL",&
   "CO(E4)+C2O4^++ANY_NEUTRAL=>C2O3^++CO2+ANY_NEUTRAL","CO(E1)+C^+=>CO^++C                               ",&
   "CO(E2)+C^+=>CO^++C                               ","CO(E3)+C^+=>CO^++C                               "/
  data reaction_sign(325:360) &
  /"CO(E4)+C^+=>CO^++C                               ","O2+CO2^+=>O2^++CO2                               ",&
   "O2+CO^+=>O2^++CO                                 ","O2+C2O2^+=>CO+CO+O2^+                            ",&
   "O2+C^+=>CO+O^+                                   ","O2+C^+=>CO^++O                                   ",&
   "O^++O2=>O2^++O                                   ","O2^++O2+ANY_NEUTRAL=>O4^++ANY_NEUTRAL            ",&
   "O2^-+O2+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL            ","O^-+O2=>O3+E                                     ",&
   "O^-+O2+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL             ","O2^-+O2=>O2+O2+E                                 ",&
   "O2+O3^-=>O2+O3+E                                 ","O2(V1)+CO2^+=>O2^++CO2                           ",&
   "O2(V2)+CO2^+=>O2^++CO2                           ","O2(V3)+CO2^+=>O2^++CO2                           ",&
   "O2(V1)+CO^+=>O2^++CO                             ","O2(V2)+CO^+=>O2^++CO                             ",&
   "O2(V3)+CO^+=>O2^++CO                             ","O2(V1)+C2O2^+=>CO+CO+O2^+                        ",&
   "O2(V2)+C2O2^+=>CO+CO+O2^+                        ","O2(V3)+C2O2^+=>CO+CO+O2^+                        ",&
   "O2(V1)+C^+=>CO+O^+                               ","O2(V2)+C^+=>CO+O^+                               ",&
   "O2(V3)+C^+=>CO+O^+                               ","O2(V1)+C^+=>CO^++O                               ",&
   "O2(V2)+C^+=>CO^++O                               ","O2(V3)+C^+=>CO^++O                               ",&
   "O2(V1)+O^+=>O2^++O                               ","O2(V2)+O^+=>O2^++O                               ",&
   "O2(V3)+O^+=>O2^++O                               ","O2(V1)+O2^++ANY_NEUTRAL=>O4^++ANY_NEUTRAL        ",&
   "O2(V2)+O2^++ANY_NEUTRAL=>O4^++ANY_NEUTRAL        ","O2(V3)+O2^++ANY_NEUTRAL=>O4^++ANY_NEUTRAL        ",&
   "O2(V1)+O2^-+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL        ","O2(V2)+O2^-+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL        "/
  data reaction_sign(361:396) &
  /"O2(V3)+O2^-+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL        ","O2(V1)+O^-=>O3+E                                 ",&
   "O2(V2)+O^-=>O3+E                                 ","O2(V3)+O^-=>O3+E                                 ",&
   "O2(V1)+O^-+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL         ","O2(V2)+O^-+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL         ",&
   "O2(V3)+O^-+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL         ","O2(V1)+O2^-=>O2(V1)+O2+E                         ",&
   "O2(V2)+O2^-=>O2(V2)+O2+E                         ","O2(V3)+O2^-=>O2(V3)+O2+E                         ",&
   "O2(V1)+O3^-=>O2(V1)+O3+E                         ","O2(V2)+O3^-=>O2(V2)+O3+E                         ",&
   "O2(V3)+O3^-=>O2(V3)+O3+E                         ","O2(E1)+CO2^+=>O2^++CO2                           ",&
   "O2(E2)+CO2^+=>O2^++CO2                           ","O2(E1)+CO^+=>O2^++CO                             ",&
   "O2(E2)+CO^+=>O2^++CO                             ","O2(E1)+C2O2^+=>CO+CO+O2^+                        ",&
   "O2(E2)+C2O2^+=>CO+CO+O2^+                        ","O2(E1)+C^+=>CO+O^+                               ",&
   "O2(E2)+C^+=>CO+O^+                               ","O2(E1)+C^+=>CO^++O                               ",&
   "O2(E2)+C^+=>CO^++O                               ","O2(E1)+O^+=>O2^++O                               ",&
   "O2(E2)+O^+=>O2^++O                               ","O2(E1)+O2^++ANY_NEUTRAL=>O4^++ANY_NEUTRAL        ",&
   "O2(E2)+O2^++ANY_NEUTRAL=>O4^++ANY_NEUTRAL        ","O2(E1)+O2^-+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL        ",&
   "O2(E2)+O2^-+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL        ","O2(E1)+O^-=>O3+E                                 ",&
   "O2(E2)+O^-=>O3+E                                 ","O2(E1)+O^-+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL         ",&
   "O2(E2)+O^-+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL         ","O2(E1)+O2^-=>O2(E1)+O2+E                         ",&
   "O2(E2)+O2^-=>O2(E2)+O2+E                         ","O2(E1)+O3^-=>O2(E1)+O3+E                         "/
  data reaction_sign(397:432) &
  /"O2(E2)+O3^-=>O2(E2)+O3+E                         ","CO^++C=>CO+C^+                                   ",&
   "O2^++C=>CO^++O                                   ","O2^++C=>C^++O2                                   ",&
   "C2^++C=>C2+C^+                                   ","O+CO2^+=>O2^++CO                                 ",&
   "O+CO2^+=>O^++CO2                                 ","CO3^-+CO2^+=>2CO2+O                              ",&
   "CO4^-+CO2^+=>2CO2+O2                             ","O2^-+CO2^+=>CO+O2+O                              ",&
   "O+CO^+=>CO+O^+                                   ","C2O2^++ANY_NEUTRAL=>CO^++CO+ANY_NEUTRAL          ",&
   "CO3^-+C2O2^+=>CO2+2CO+O                          ","CO4^-+C2O2^+=>CO2+2CO+O2                         ",&
   "O2^-+C2O2^+=>2CO+O2                              ","CO3^-+C2O3^+=>2CO2+CO+O                          ",&
   "CO4^-+C2O3^+=>2CO2+CO+O2                         ","O2^-+C2O3^+=>CO2+CO+O2                           ",&
   "C2O4^++ANY_NEUTRAL=>CO2^++CO2+ANY_NEUTRAL        ","CO3^-+C2O4^+=>3CO2+O                             ",&
   "CO4^-+C2O4^+=>3CO2+O2                            ","O2^-+C2O4^+=>2CO2+O2                             ",&
   "O2^++CO3^-=>CO2+O2+O                             ","O+CO3^-=>CO2+O2^-                                ",&
   "O2^++CO4^-=>CO2+O2+O2                            ","O+CO4^-=>CO3^-+O2                                ",&
   "O+CO4^-=>CO2+O2+O^-                              ","O+CO4^-=>CO2+O3^-                                ",&
   "O3+CO4^-=>CO2+O3^-+O2                            ","O^-+O3=>O3^-+O                                   ",&
   "O^-+O3=>O2+O2+E                                  ","O2^-+O3=>O3^-+O2                                 ",&
   "O3^-+O3=>O2+O2+O2+E                              ","O^++O3=>O2^++O2                                  ",&
   "O^++O+ANY_NEUTRAL=>O2^++ANY_NEUTRAL              ","O^-+O=>O2+E                                      "/
  data reaction_sign(433:468) &
  /"O2^-+O=>O2+O^-                                   ","O2^-+O=>O3+E                                     ",&
   "O3^-+O=>O3+O^-                                   ","O3^-+O=>O2+O2+E                                  ",&
   "O3^-+O=>O2^-+O2                                  ","O4^-+O=>O3^-+O2                                  ",&
   "O4^-+O=>O^-+O2+O2                                ","O4^-+O=>O2^-+O3                                  ",&
   "O^++O2^-+ANY_NEUTRAL=>O3+ANY_NEUTRAL             ","O^++O2^-=>O+O2                                   ",&
   "O2^++O2^-=>O2+O2                                 ","O2^++O2^-=>O2+O+O                                ",&
   "O2^++O2^-+ANY_NEUTRAL=>O2+O2+ANY_NEUTRAL         ","O2^-+ANY_NEUTRAL=>O2+ANY_NEUTRAL+E               ",&
   "O2^++O3^-=>O2+O3                                 ","O2^++O3^-=>O+O+O3                                ",&
   "O^++O3^-=>O3+O                                   ","O3^-+ANY_NEUTRAL=>O3+E+ANY_NEUTRAL               ",&
   "O^++O^-=>O+O                                     ","O^++O^-+ANY_NEUTRAL=>O2+ANY_NEUTRAL              ",&
   "O2^++O^-=>O2+O                                   ","O2^++O^-=>O+O+O                                  ",&
   "O2^++O^-+ANY_NEUTRAL=>O3+ANY_NEUTRAL             ","ANY_NEUTRAL+O^-=>O+ANY_NEUTRAL+E                 ",&
   "ANY_NEUTRAL+O4^-=>O2^-+O2+ANY_NEUTRAL            ","O4^++ANY_NEUTRAL=>O2^++O2+ANY_NEUTRAL            ",&
   "CO2(VA)+CO2=>CO2+CO2                             ","CO2(VA)+CO=>CO2+CO                               ",&
   "CO2(VA)+O2=>CO2+O2                               ","CO2(VB)+CO2=>CO2+CO2                             ",&
   "CO2(VB)+CO=>CO2+CO                               ","CO2(VB)+O2=>CO2+O2                               ",&
   "CO2(VB)+CO2=>CO2(VA)+CO2                         ","CO2(VB)+CO=>CO2(VA)+CO                           ",&
   "CO2(VB)+O2=>CO2(VA)+O2                           ","CO2(V1)+CO2=>CO2(VB)+CO2                         "/
  data reaction_sign(469:503) &
  /"CO2(V1)+CO=>CO2(VB)+CO                           ","CO2(V1)+O2=>CO2(VB)+O2                           ",&
   "CO2(V1)+CO2=>CO2(VD)+CO2                         ","CO2(V1)+CO=>CO2(VD)+CO                           ",&
   "CO2(V1)+O2=>CO2(VD)+O2                           ","CO2(V1)+CO2=>CO2(VA)+CO2(VB)                     ",&
   "CO2(V1)+CO2=>CO2(VA)+CO2                         ","CO2(V1)+CO=>CO2(VA)+CO                           ",&
   "CO2(V1)+O2=>CO2(VA)+O2                           ","CO2(VD)+CO2=>CO2(VB)+CO2                         ",&
   "CO2(VD)+CO=>CO2(VB)+CO                           ","CO2(VD)+O2=>CO2(VB)+O2                           ",&
   "CO2(VD)+CO2=>CO2(VA)+CO2                         ","CO2(VD)+CO=>CO2(VA)+CO                           ",&
   "CO2(VD)+O2=>CO2(VA)+O2                           ","CO(V1)+CO2=>CO+CO2                               ",&
   "CO(V1)+CO=>CO+CO                                 ","CO(V1)+O2=>CO+O2                                 ",&
   "O2(V1)+CO2=>O2+CO2                               ","O2(V1)+CO=>O2+CO                                 ",&
   "O2(V1)+O2=>O2+O2                                 ","O2(V2)+CO2=>O2+CO2                               ",&
   "O2(V2)+CO=>O2+CO                                 ","O2(V2)+O2=>O2+O2                                 ",&
   "O2(V3)+CO2=>O2+CO2                               ","O2(V3)+CO=>O2+CO                                 ",&
   "O2(V3)+O2=>O2+O2                                 ","CO2(E1)+ANY_NEUTRAL=>CO2+ANY_NEUTRAL             ",&
   "CO2(E2)+ANY_NEUTRAL=>CO2+ANY_NEUTRAL             ","CO(E1)+ANY_NEUTRAL=>CO+ANY_NEUTRAL               ",&
   "CO(E2)+ANY_NEUTRAL=>CO+ANY_NEUTRAL               ","CO(E3)+ANY_NEUTRAL=>CO+ANY_NEUTRAL               ",&
   "CO(E4)+ANY_NEUTRAL=>CO+ANY_NEUTRAL               ","O2(E1)+ANY_NEUTRAL=>O2+ANY_NEUTRAL               ",&
   "O2(E2)+ANY_NEUTRAL=>O2+ANY_NEUTRAL               "/
  data bolsig_species(1:bolsig_species_max) &
  /"CO2    ","CO2(VA)","CO2(VB)","CO2(VD)","CO2(V1)","CO2(E1)","CO2(E2)","CO     ","CO(V1) ","CO(E1) ","CO(E2) ","CO(E3) ",&
   "CO(E4) ","O2     ","O2(V1) ","O2(V2) ","O2(V3) ","O2(E1) ","O2(E2) ","O3     ","O      ","C      ","C2     "/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A7)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A49)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(541x,42(1x,i9))
312 format(A3,1x,A49,1x,42(1x,A9))
313 format(i3,1x,A49,1x,42(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(503(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,42(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,503(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,42(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(504(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(01,001) = + reac_rate_local(001) 
  reac_source_local(02,001) = - reac_rate_local(001) 
  reac_source_local(26,001) = + reac_rate_local(001) 
  reac_source_local(01,002) = + reac_rate_local(002) 
  reac_source_local(02,002) = - reac_rate_local(002) 
  reac_source_local(25,002) = + reac_rate_local(002) 
  reac_source_local(28,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = + reac_rate_local(003) 
  reac_source_local(02,003) = - reac_rate_local(003) 
  reac_source_local(15,003) = + reac_rate_local(003) 
  reac_source_local(32,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = + reac_rate_local(004) 
  reac_source_local(02,004) = - reac_rate_local(004) 
  reac_source_local(09,004) = + reac_rate_local(004) 
  reac_source_local(37,004) = + reac_rate_local(004) 
  reac_source_local(01,005) = + reac_rate_local(005) 
  reac_source_local(02,005) = - reac_rate_local(005) 
  reac_source_local(22,005) = + reac_rate_local(005) 
  reac_source_local(36,005) = + reac_rate_local(005) 
  reac_source_local(01,006) = - reac_rate_local(006) 
  reac_source_local(02,006) = - reac_rate_local(006) 
  reac_source_local(09,006) = + reac_rate_local(006) 
  reac_source_local(39,006) = + reac_rate_local(006) 
  reac_source_local(02,007) = - reac_rate_local(007) 
  reac_source_local(09,007) = + reac_rate_local(007) 
  reac_source_local(25,007) = + reac_rate_local(007) 
  reac_source_local(02,008) = - reac_rate_local(008) 
  reac_source_local(03,008) = + reac_rate_local(008) 
  reac_source_local(02,009) = - reac_rate_local(009) 
  reac_source_local(04,009) = + reac_rate_local(009) 
  reac_source_local(02,010) = - reac_rate_local(010) 
  reac_source_local(05,010) = + reac_rate_local(010) 
  reac_source_local(02,011) = - reac_rate_local(011) 
  reac_source_local(06,011) = + reac_rate_local(011) 
  reac_source_local(02,012) = - reac_rate_local(012) 
  reac_source_local(07,012) = + reac_rate_local(012) 
  reac_source_local(02,013) = - reac_rate_local(013) 
  reac_source_local(08,013) = + reac_rate_local(013) 
  reac_source_local(02,014) = + reac_rate_local(014) 
  reac_source_local(07,014) = - reac_rate_local(014) 
  reac_source_local(02,015) = + reac_rate_local(015) 
  reac_source_local(08,015) = - reac_rate_local(015) 
  reac_source_local(03,016) = - reac_rate_local(016) 
  reac_source_local(09,016) = + reac_rate_local(016) 
  reac_source_local(25,016) = + reac_rate_local(016) 
  reac_source_local(04,017) = - reac_rate_local(017) 
  reac_source_local(09,017) = + reac_rate_local(017) 
  reac_source_local(25,017) = + reac_rate_local(017) 
  reac_source_local(06,018) = - reac_rate_local(018) 
  reac_source_local(09,018) = + reac_rate_local(018) 
  reac_source_local(25,018) = + reac_rate_local(018) 
  reac_source_local(05,019) = - reac_rate_local(019) 
  reac_source_local(09,019) = + reac_rate_local(019) 
  reac_source_local(25,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(03,020) = - reac_rate_local(020) 
  reac_source_local(26,020) = + reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(04,021) = - reac_rate_local(021) 
  reac_source_local(26,021) = + reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(06,022) = - reac_rate_local(022) 
  reac_source_local(26,022) = + reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(05,023) = - reac_rate_local(023) 
  reac_source_local(26,023) = + reac_rate_local(023) 
  reac_source_local(01,024) = + reac_rate_local(024) 
  reac_source_local(03,024) = - reac_rate_local(024) 
  reac_source_local(25,024) = + reac_rate_local(024) 
  reac_source_local(28,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(04,025) = - reac_rate_local(025) 
  reac_source_local(25,025) = + reac_rate_local(025) 
  reac_source_local(28,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(06,026) = - reac_rate_local(026) 
  reac_source_local(25,026) = + reac_rate_local(026) 
  reac_source_local(28,026) = + reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(05,027) = - reac_rate_local(027) 
  reac_source_local(25,027) = + reac_rate_local(027) 
  reac_source_local(28,027) = + reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(03,028) = - reac_rate_local(028) 
  reac_source_local(15,028) = + reac_rate_local(028) 
  reac_source_local(32,028) = + reac_rate_local(028) 
  reac_source_local(01,029) = + reac_rate_local(029) 
  reac_source_local(04,029) = - reac_rate_local(029) 
  reac_source_local(15,029) = + reac_rate_local(029) 
  reac_source_local(32,029) = + reac_rate_local(029) 
  reac_source_local(01,030) = + reac_rate_local(030) 
  reac_source_local(06,030) = - reac_rate_local(030) 
  reac_source_local(15,030) = + reac_rate_local(030) 
  reac_source_local(32,030) = + reac_rate_local(030) 
  reac_source_local(01,031) = + reac_rate_local(031) 
  reac_source_local(05,031) = - reac_rate_local(031) 
  reac_source_local(15,031) = + reac_rate_local(031) 
  reac_source_local(32,031) = + reac_rate_local(031) 
  reac_source_local(01,032) = + reac_rate_local(032) 
  reac_source_local(03,032) = - reac_rate_local(032) 
  reac_source_local(09,032) = + reac_rate_local(032) 
  reac_source_local(37,032) = + reac_rate_local(032) 
  reac_source_local(01,033) = + reac_rate_local(033) 
  reac_source_local(04,033) = - reac_rate_local(033) 
  reac_source_local(09,033) = + reac_rate_local(033) 
  reac_source_local(37,033) = + reac_rate_local(033) 
  reac_source_local(01,034) = + reac_rate_local(034) 
  reac_source_local(06,034) = - reac_rate_local(034) 
  reac_source_local(09,034) = + reac_rate_local(034) 
  reac_source_local(37,034) = + reac_rate_local(034) 
  reac_source_local(01,035) = + reac_rate_local(035) 
  reac_source_local(05,035) = - reac_rate_local(035) 
  reac_source_local(09,035) = + reac_rate_local(035) 
  reac_source_local(37,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = - reac_rate_local(036) 
  reac_source_local(03,036) = - reac_rate_local(036) 
  reac_source_local(09,036) = + reac_rate_local(036) 
  reac_source_local(39,036) = + reac_rate_local(036) 
  reac_source_local(01,037) = - reac_rate_local(037) 
  reac_source_local(04,037) = - reac_rate_local(037) 
  reac_source_local(09,037) = + reac_rate_local(037) 
  reac_source_local(39,037) = + reac_rate_local(037) 
  reac_source_local(01,038) = - reac_rate_local(038) 
  reac_source_local(06,038) = - reac_rate_local(038) 
  reac_source_local(09,038) = + reac_rate_local(038) 
  reac_source_local(39,038) = + reac_rate_local(038) 
  reac_source_local(01,039) = - reac_rate_local(039) 
  reac_source_local(05,039) = - reac_rate_local(039) 
  reac_source_local(09,039) = + reac_rate_local(039) 
  reac_source_local(39,039) = + reac_rate_local(039) 
  reac_source_local(03,040) = - reac_rate_local(040) 
  reac_source_local(07,040) = + reac_rate_local(040) 
  reac_source_local(04,041) = - reac_rate_local(041) 
  reac_source_local(07,041) = + reac_rate_local(041) 
  reac_source_local(06,042) = - reac_rate_local(042) 
  reac_source_local(07,042) = + reac_rate_local(042) 
  reac_source_local(05,043) = - reac_rate_local(043) 
  reac_source_local(07,043) = + reac_rate_local(043) 
  reac_source_local(03,044) = - reac_rate_local(044) 
  reac_source_local(08,044) = + reac_rate_local(044) 
  reac_source_local(04,045) = - reac_rate_local(045) 
  reac_source_local(08,045) = + reac_rate_local(045) 
  reac_source_local(06,046) = - reac_rate_local(046) 
  reac_source_local(08,046) = + reac_rate_local(046) 
  reac_source_local(05,047) = - reac_rate_local(047) 
  reac_source_local(08,047) = + reac_rate_local(047) 
  reac_source_local(01,048) = + reac_rate_local(048) 
  reac_source_local(03,048) = - reac_rate_local(048) 
  reac_source_local(22,048) = + reac_rate_local(048) 
  reac_source_local(36,048) = + reac_rate_local(048) 
  reac_source_local(01,049) = + reac_rate_local(049) 
  reac_source_local(04,049) = - reac_rate_local(049) 
  reac_source_local(22,049) = + reac_rate_local(049) 
  reac_source_local(36,049) = + reac_rate_local(049) 
  reac_source_local(01,050) = + reac_rate_local(050) 
  reac_source_local(06,050) = - reac_rate_local(050) 
  reac_source_local(22,050) = + reac_rate_local(050) 
  reac_source_local(36,050) = + reac_rate_local(050) 
  reac_source_local(01,051) = + reac_rate_local(051) 
  reac_source_local(05,051) = - reac_rate_local(051) 
  reac_source_local(22,051) = + reac_rate_local(051) 
  reac_source_local(36,051) = + reac_rate_local(051) 
  reac_source_local(02,052) = + reac_rate_local(052) 
  reac_source_local(03,052) = - reac_rate_local(052) 
  reac_source_local(02,053) = + reac_rate_local(053) 
  reac_source_local(04,053) = - reac_rate_local(053) 
  reac_source_local(02,054) = + reac_rate_local(054) 
  reac_source_local(06,054) = - reac_rate_local(054) 
  reac_source_local(02,055) = + reac_rate_local(055) 
  reac_source_local(05,055) = - reac_rate_local(055) 
  reac_source_local(07,056) = - reac_rate_local(056) 
  reac_source_local(09,056) = + reac_rate_local(056) 
  reac_source_local(25,056) = + reac_rate_local(056) 
  reac_source_local(08,057) = - reac_rate_local(057) 
  reac_source_local(09,057) = + reac_rate_local(057) 
  reac_source_local(25,057) = + reac_rate_local(057) 
  reac_source_local(01,058) = + reac_rate_local(058) 
  reac_source_local(07,058) = - reac_rate_local(058) 
  reac_source_local(26,058) = + reac_rate_local(058) 
  reac_source_local(01,059) = + reac_rate_local(059) 
  reac_source_local(08,059) = - reac_rate_local(059) 
  reac_source_local(26,059) = + reac_rate_local(059) 
  reac_source_local(01,060) = + reac_rate_local(060) 
  reac_source_local(07,060) = - reac_rate_local(060) 
  reac_source_local(25,060) = + reac_rate_local(060) 
  reac_source_local(28,060) = + reac_rate_local(060) 
  reac_source_local(01,061) = + reac_rate_local(061) 
  reac_source_local(08,061) = - reac_rate_local(061) 
  reac_source_local(25,061) = + reac_rate_local(061) 
  reac_source_local(28,061) = + reac_rate_local(061) 
  reac_source_local(01,062) = + reac_rate_local(062) 
  reac_source_local(07,062) = - reac_rate_local(062) 
  reac_source_local(15,062) = + reac_rate_local(062) 
  reac_source_local(32,062) = + reac_rate_local(062) 
  reac_source_local(01,063) = + reac_rate_local(063) 
  reac_source_local(08,063) = - reac_rate_local(063) 
  reac_source_local(15,063) = + reac_rate_local(063) 
  reac_source_local(32,063) = + reac_rate_local(063) 
  reac_source_local(01,064) = + reac_rate_local(064) 
  reac_source_local(07,064) = - reac_rate_local(064) 
  reac_source_local(09,064) = + reac_rate_local(064) 
  reac_source_local(37,064) = + reac_rate_local(064) 
  reac_source_local(01,065) = + reac_rate_local(065) 
  reac_source_local(08,065) = - reac_rate_local(065) 
  reac_source_local(09,065) = + reac_rate_local(065) 
  reac_source_local(37,065) = + reac_rate_local(065) 
  reac_source_local(01,066) = + reac_rate_local(066) 
  reac_source_local(07,066) = - reac_rate_local(066) 
  reac_source_local(22,066) = + reac_rate_local(066) 
  reac_source_local(36,066) = + reac_rate_local(066) 
  reac_source_local(01,067) = + reac_rate_local(067) 
  reac_source_local(08,067) = - reac_rate_local(067) 
  reac_source_local(22,067) = + reac_rate_local(067) 
  reac_source_local(36,067) = + reac_rate_local(067) 
  reac_source_local(01,068) = + reac_rate_local(068) 
  reac_source_local(09,068) = - reac_rate_local(068) 
  reac_source_local(28,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(09,069) = - reac_rate_local(069) 
  reac_source_local(25,069) = + reac_rate_local(069) 
  reac_source_local(32,069) = + reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(09,070) = - reac_rate_local(070) 
  reac_source_local(22,070) = + reac_rate_local(070) 
  reac_source_local(37,070) = + reac_rate_local(070) 
  reac_source_local(01,071) = - reac_rate_local(071) 
  reac_source_local(09,071) = - reac_rate_local(071) 
  reac_source_local(22,071) = + reac_rate_local(071) 
  reac_source_local(39,071) = + reac_rate_local(071) 
  reac_source_local(09,072) = - reac_rate_local(072) 
  reac_source_local(22,072) = + reac_rate_local(072) 
  reac_source_local(25,072) = + reac_rate_local(072) 
  reac_source_local(09,073) = - reac_rate_local(073) 
  reac_source_local(10,073) = + reac_rate_local(073) 
  reac_source_local(09,074) = - reac_rate_local(074) 
  reac_source_local(11,074) = + reac_rate_local(074) 
  reac_source_local(09,075) = - reac_rate_local(075) 
  reac_source_local(12,075) = + reac_rate_local(075) 
  reac_source_local(09,076) = - reac_rate_local(076) 
  reac_source_local(13,076) = + reac_rate_local(076) 
  reac_source_local(09,077) = - reac_rate_local(077) 
  reac_source_local(14,077) = + reac_rate_local(077) 
  reac_source_local(09,078) = + reac_rate_local(078) 
  reac_source_local(11,078) = - reac_rate_local(078) 
  reac_source_local(09,079) = + reac_rate_local(079) 
  reac_source_local(12,079) = - reac_rate_local(079) 
  reac_source_local(09,080) = + reac_rate_local(080) 
  reac_source_local(13,080) = - reac_rate_local(080) 
  reac_source_local(09,081) = + reac_rate_local(081) 
  reac_source_local(14,081) = - reac_rate_local(081) 
  reac_source_local(10,082) = - reac_rate_local(082) 
  reac_source_local(22,082) = + reac_rate_local(082) 
  reac_source_local(25,082) = + reac_rate_local(082) 
  reac_source_local(01,083) = + reac_rate_local(083) 
  reac_source_local(10,083) = - reac_rate_local(083) 
  reac_source_local(28,083) = + reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(10,084) = - reac_rate_local(084) 
  reac_source_local(25,084) = + reac_rate_local(084) 
  reac_source_local(32,084) = + reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(10,085) = - reac_rate_local(085) 
  reac_source_local(22,085) = + reac_rate_local(085) 
  reac_source_local(37,085) = + reac_rate_local(085) 
  reac_source_local(01,086) = - reac_rate_local(086) 
  reac_source_local(10,086) = - reac_rate_local(086) 
  reac_source_local(22,086) = + reac_rate_local(086) 
  reac_source_local(39,086) = + reac_rate_local(086) 
  reac_source_local(10,087) = - reac_rate_local(087) 
  reac_source_local(11,087) = + reac_rate_local(087) 
  reac_source_local(10,088) = - reac_rate_local(088) 
  reac_source_local(12,088) = + reac_rate_local(088) 
  reac_source_local(10,089) = - reac_rate_local(089) 
  reac_source_local(13,089) = + reac_rate_local(089) 
  reac_source_local(10,090) = - reac_rate_local(090) 
  reac_source_local(14,090) = + reac_rate_local(090) 
  reac_source_local(09,091) = + reac_rate_local(091) 
  reac_source_local(10,091) = - reac_rate_local(091) 
  reac_source_local(11,092) = - reac_rate_local(092) 
  reac_source_local(22,092) = + reac_rate_local(092) 
  reac_source_local(25,092) = + reac_rate_local(092) 
  reac_source_local(12,093) = - reac_rate_local(093) 
  reac_source_local(22,093) = + reac_rate_local(093) 
  reac_source_local(25,093) = + reac_rate_local(093) 
  reac_source_local(13,094) = - reac_rate_local(094) 
  reac_source_local(22,094) = + reac_rate_local(094) 
  reac_source_local(25,094) = + reac_rate_local(094) 
  reac_source_local(14,095) = - reac_rate_local(095) 
  reac_source_local(22,095) = + reac_rate_local(095) 
  reac_source_local(25,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(11,096) = - reac_rate_local(096) 
  reac_source_local(28,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(12,097) = - reac_rate_local(097) 
  reac_source_local(28,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(13,098) = - reac_rate_local(098) 
  reac_source_local(28,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(14,099) = - reac_rate_local(099) 
  reac_source_local(28,099) = + reac_rate_local(099) 
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(11,100) = - reac_rate_local(100) 
  reac_source_local(25,100) = + reac_rate_local(100) 
  reac_source_local(32,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(12,101) = - reac_rate_local(101) 
  reac_source_local(25,101) = + reac_rate_local(101) 
  reac_source_local(32,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(13,102) = - reac_rate_local(102) 
  reac_source_local(25,102) = + reac_rate_local(102) 
  reac_source_local(32,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(14,103) = - reac_rate_local(103) 
  reac_source_local(25,103) = + reac_rate_local(103) 
  reac_source_local(32,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(11,104) = - reac_rate_local(104) 
  reac_source_local(22,104) = + reac_rate_local(104) 
  reac_source_local(37,104) = + reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(12,105) = - reac_rate_local(105) 
  reac_source_local(22,105) = + reac_rate_local(105) 
  reac_source_local(37,105) = + reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(13,106) = - reac_rate_local(106) 
  reac_source_local(22,106) = + reac_rate_local(106) 
  reac_source_local(37,106) = + reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(14,107) = - reac_rate_local(107) 
  reac_source_local(22,107) = + reac_rate_local(107) 
  reac_source_local(37,107) = + reac_rate_local(107) 
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(22,108) = - reac_rate_local(108) 
  reac_source_local(32,108) = + reac_rate_local(108) 
  reac_source_local(22,109) = + reac_rate_local(109) * 2.d0
  reac_source_local(23,109) = - reac_rate_local(109) 
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(23,110) = - reac_rate_local(110) 
  reac_source_local(33,110) = + reac_rate_local(110) 
  reac_source_local(15,111) = - reac_rate_local(111) 
  reac_source_local(25,111) = + reac_rate_local(111) * 2.d0
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(15,112) = - reac_rate_local(112) 
  reac_source_local(36,112) = + reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(15,113) = - reac_rate_local(113) 
  reac_source_local(25,113) = + reac_rate_local(113) 
  reac_source_local(37,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = - reac_rate_local(114) 
  reac_source_local(15,114) = - reac_rate_local(114) 
  reac_source_local(25,114) = + reac_rate_local(114) 
  reac_source_local(39,114) = + reac_rate_local(114) 
  reac_source_local(15,115) = - reac_rate_local(115) 
  reac_source_local(16,115) = + reac_rate_local(115) 
  reac_source_local(15,116) = - reac_rate_local(116) 
  reac_source_local(17,116) = + reac_rate_local(116) 
  reac_source_local(15,117) = - reac_rate_local(117) 
  reac_source_local(18,117) = + reac_rate_local(117) 
  reac_source_local(15,118) = - reac_rate_local(118) 
  reac_source_local(19,118) = + reac_rate_local(118) 
  reac_source_local(15,119) = - reac_rate_local(119) 
  reac_source_local(20,119) = + reac_rate_local(119) 
  reac_source_local(15,120) = + reac_rate_local(120) 
  reac_source_local(19,120) = - reac_rate_local(120) 
  reac_source_local(15,121) = + reac_rate_local(121) 
  reac_source_local(20,121) = - reac_rate_local(121) 
  reac_source_local(15,122) = + reac_rate_local(122) 
  reac_source_local(19,122) = - reac_rate_local(122) 
  reac_source_local(15,123) = + reac_rate_local(123) 
  reac_source_local(20,123) = - reac_rate_local(123) 
  reac_source_local(16,124) = - reac_rate_local(124) 
  reac_source_local(19,124) = + reac_rate_local(124) 
  reac_source_local(17,125) = - reac_rate_local(125) 
  reac_source_local(19,125) = + reac_rate_local(125) 
  reac_source_local(18,126) = - reac_rate_local(126) 
  reac_source_local(19,126) = + reac_rate_local(126) 
  reac_source_local(16,127) = - reac_rate_local(127) 
  reac_source_local(20,127) = + reac_rate_local(127) 
  reac_source_local(17,128) = - reac_rate_local(128) 
  reac_source_local(20,128) = + reac_rate_local(128) 
  reac_source_local(18,129) = - reac_rate_local(129) 
  reac_source_local(20,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(16,130) = - reac_rate_local(130) 
  reac_source_local(36,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(17,131) = - reac_rate_local(131) 
  reac_source_local(36,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(18,132) = - reac_rate_local(132) 
  reac_source_local(36,132) = + reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(16,133) = - reac_rate_local(133) 
  reac_source_local(25,133) = + reac_rate_local(133) 
  reac_source_local(37,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(17,134) = - reac_rate_local(134) 
  reac_source_local(25,134) = + reac_rate_local(134) 
  reac_source_local(37,134) = + reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(18,135) = - reac_rate_local(135) 
  reac_source_local(25,135) = + reac_rate_local(135) 
  reac_source_local(37,135) = + reac_rate_local(135) 
  reac_source_local(16,136) = - reac_rate_local(136) 
  reac_source_local(25,136) = + reac_rate_local(136) * 2.d0
  reac_source_local(17,137) = - reac_rate_local(137) 
  reac_source_local(25,137) = + reac_rate_local(137) * 2.d0
  reac_source_local(18,138) = - reac_rate_local(138) 
  reac_source_local(25,138) = + reac_rate_local(138) * 2.d0
  reac_source_local(01,139) = - reac_rate_local(139) 
  reac_source_local(16,139) = - reac_rate_local(139) 
  reac_source_local(25,139) = + reac_rate_local(139) 
  reac_source_local(39,139) = + reac_rate_local(139) 
  reac_source_local(01,140) = - reac_rate_local(140) 
  reac_source_local(17,140) = - reac_rate_local(140) 
  reac_source_local(25,140) = + reac_rate_local(140) 
  reac_source_local(39,140) = + reac_rate_local(140) 
  reac_source_local(01,141) = - reac_rate_local(141) 
  reac_source_local(18,141) = - reac_rate_local(141) 
  reac_source_local(25,141) = + reac_rate_local(141) 
  reac_source_local(39,141) = + reac_rate_local(141) 
  reac_source_local(15,142) = + reac_rate_local(142) 
  reac_source_local(16,142) = - reac_rate_local(142) 
  reac_source_local(15,143) = + reac_rate_local(143) 
  reac_source_local(17,143) = - reac_rate_local(143) 
  reac_source_local(15,144) = + reac_rate_local(144) 
  reac_source_local(18,144) = - reac_rate_local(144) 
  reac_source_local(19,145) = - reac_rate_local(145) 
  reac_source_local(25,145) = + reac_rate_local(145) * 2.d0
  reac_source_local(20,146) = - reac_rate_local(146) 
  reac_source_local(25,146) = + reac_rate_local(146) * 2.d0
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(19,147) = - reac_rate_local(147) 
  reac_source_local(36,147) = + reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(20,148) = - reac_rate_local(148) 
  reac_source_local(36,148) = + reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(19,149) = - reac_rate_local(149) 
  reac_source_local(25,149) = + reac_rate_local(149) 
  reac_source_local(37,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(20,150) = - reac_rate_local(150) 
  reac_source_local(25,150) = + reac_rate_local(150) 
  reac_source_local(37,150) = + reac_rate_local(150) 
  reac_source_local(15,151) = + reac_rate_local(151) 
  reac_source_local(24,151) = - reac_rate_local(151) 
  reac_source_local(25,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(24,152) = - reac_rate_local(152) 
  reac_source_local(25,152) = + reac_rate_local(152) 
  reac_source_local(36,152) = + reac_rate_local(152) 
  reac_source_local(24,153) = - reac_rate_local(153) 
  reac_source_local(25,153) = + reac_rate_local(153) 
  reac_source_local(37,153) = + reac_rate_local(153) 
  reac_source_local(39,153) = + reac_rate_local(153) 
  reac_source_local(01,154) = - reac_rate_local(154) 
  reac_source_local(15,154) = + reac_rate_local(154) 
  reac_source_local(24,154) = - reac_rate_local(154) 
  reac_source_local(39,154) = + reac_rate_local(154) 
  reac_source_local(01,155) = - reac_rate_local(155) 
  reac_source_local(24,155) = - reac_rate_local(155) 
  reac_source_local(25,155) = + reac_rate_local(155) 
  reac_source_local(40,155) = + reac_rate_local(155) 
  reac_source_local(01,156) = + reac_rate_local(156) 
  reac_source_local(25,156) = - reac_rate_local(156) 
  reac_source_local(37,156) = + reac_rate_local(156) 
  reac_source_local(01,157) = - reac_rate_local(157) 
  reac_source_local(09,157) = + reac_rate_local(157) 
  reac_source_local(25,157) = + reac_rate_local(157) 
  reac_source_local(26,157) = - reac_rate_local(157) 
  reac_source_local(01,158) = - reac_rate_local(158) 
  reac_source_local(15,158) = + reac_rate_local(158) 
  reac_source_local(22,158) = + reac_rate_local(158) 
  reac_source_local(26,158) = - reac_rate_local(158) 
  reac_source_local(01,159) = - reac_rate_local(159) 
  reac_source_local(02,159) = + reac_rate_local(159) 
  reac_source_local(15,159) = + reac_rate_local(159) 
  reac_source_local(27,159) = - reac_rate_local(159) 
  reac_source_local(01,160) = - reac_rate_local(160) 
  reac_source_local(22,160) = + reac_rate_local(160) 
  reac_source_local(25,160) = + reac_rate_local(160) 
  reac_source_local(28,160) = - reac_rate_local(160) 
  reac_source_local(01,161) = - reac_rate_local(161) 
  reac_source_local(09,161) = + reac_rate_local(161) * 2.d0
  reac_source_local(29,161) = - reac_rate_local(161) 
  reac_source_local(01,162) = - reac_rate_local(162) 
  reac_source_local(02,162) = + reac_rate_local(162) 
  reac_source_local(09,162) = + reac_rate_local(162) 
  reac_source_local(30,162) = - reac_rate_local(162) 
  reac_source_local(01,163) = - reac_rate_local(163) 
  reac_source_local(02,163) = + reac_rate_local(163) * 2.d0
  reac_source_local(31,163) = - reac_rate_local(163) 
  reac_source_local(01,164) = - reac_rate_local(164) 
  reac_source_local(22,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(33,164) = - reac_rate_local(164) 
  reac_source_local(01,165) = - reac_rate_local(165) 
  reac_source_local(15,165) = - reac_rate_local(165) 
  reac_source_local(40,165) = + reac_rate_local(165) 
  reac_source_local(01,166) = - reac_rate_local(166) 
  reac_source_local(24,166) = - reac_rate_local(166) 
  reac_source_local(41,166) = + reac_rate_local(166) 
  reac_source_local(01,167) = - reac_rate_local(167) 
  reac_source_local(25,167) = - reac_rate_local(167) 
  reac_source_local(39,167) = + reac_rate_local(167) 
  reac_source_local(01,168) = - reac_rate_local(168) 
  reac_source_local(15,168) = + reac_rate_local(168) 
  reac_source_local(36,168) = - reac_rate_local(168) 
  reac_source_local(01,169) = - reac_rate_local(169) 
  reac_source_local(25,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(36,169) = - reac_rate_local(169) 
  reac_source_local(01,170) = - reac_rate_local(170) 
  reac_source_local(25,170) = + reac_rate_local(170) 
  reac_source_local(37,170) = - reac_rate_local(170) 
  reac_source_local(01,171) = - reac_rate_local(171) 
  reac_source_local(15,171) = + reac_rate_local(171) * 2.d0
  reac_source_local(38,171) = - reac_rate_local(171) 
  reac_source_local(02,172) = - reac_rate_local(172) 
  reac_source_local(09,172) = + reac_rate_local(172) 
  reac_source_local(15,172) = + reac_rate_local(172) 
  reac_source_local(25,172) = - reac_rate_local(172) 
  reac_source_local(02,173) = - reac_rate_local(173) 
  reac_source_local(09,173) = + reac_rate_local(173) * 2.d0
  reac_source_local(22,173) = - reac_rate_local(173) 
  reac_source_local(03,174) = - reac_rate_local(174) 
  reac_source_local(09,174) = + reac_rate_local(174) 
  reac_source_local(15,174) = + reac_rate_local(174) 
  reac_source_local(25,174) = - reac_rate_local(174) 
  reac_source_local(04,175) = - reac_rate_local(175) 
  reac_source_local(09,175) = + reac_rate_local(175) 
  reac_source_local(15,175) = + reac_rate_local(175) 
  reac_source_local(25,175) = - reac_rate_local(175) 
  reac_source_local(05,176) = - reac_rate_local(176) 
  reac_source_local(09,176) = + reac_rate_local(176) 
  reac_source_local(15,176) = + reac_rate_local(176) 
  reac_source_local(25,176) = - reac_rate_local(176) 
  reac_source_local(06,177) = - reac_rate_local(177) 
  reac_source_local(09,177) = + reac_rate_local(177) 
  reac_source_local(15,177) = + reac_rate_local(177) 
  reac_source_local(25,177) = - reac_rate_local(177) 
  reac_source_local(03,178) = - reac_rate_local(178) 
  reac_source_local(09,178) = + reac_rate_local(178) * 2.d0
  reac_source_local(22,178) = - reac_rate_local(178) 
  reac_source_local(04,179) = - reac_rate_local(179) 
  reac_source_local(09,179) = + reac_rate_local(179) * 2.d0
  reac_source_local(22,179) = - reac_rate_local(179) 
  reac_source_local(05,180) = - reac_rate_local(180) 
  reac_source_local(09,180) = + reac_rate_local(180) * 2.d0
  reac_source_local(22,180) = - reac_rate_local(180) 
  reac_source_local(06,181) = - reac_rate_local(181) 
  reac_source_local(09,181) = + reac_rate_local(181) * 2.d0
  reac_source_local(22,181) = - reac_rate_local(181) 
  reac_source_local(02,182) = + reac_rate_local(182) 
  reac_source_local(09,182) = - reac_rate_local(182) 
  reac_source_local(25,182) = - reac_rate_local(182) 
  reac_source_local(02,183) = + reac_rate_local(183) 
  reac_source_local(09,183) = - reac_rate_local(183) 
  reac_source_local(15,183) = + reac_rate_local(183) 
  reac_source_local(24,183) = - reac_rate_local(183) 
  reac_source_local(09,184) = - reac_rate_local(184) 
  reac_source_local(21,184) = + reac_rate_local(184) 
  reac_source_local(22,184) = - reac_rate_local(184) 
  reac_source_local(02,185) = + reac_rate_local(185) 
  reac_source_local(09,185) = - reac_rate_local(185) 
  reac_source_local(15,185) = - reac_rate_local(185) 
  reac_source_local(25,185) = + reac_rate_local(185) 
  reac_source_local(02,186) = + reac_rate_local(186) 
  reac_source_local(10,186) = - reac_rate_local(186) 
  reac_source_local(25,186) = - reac_rate_local(186) 
  reac_source_local(02,187) = + reac_rate_local(187) 
  reac_source_local(09,187) = - reac_rate_local(187) 
  reac_source_local(16,187) = - reac_rate_local(187) 
  reac_source_local(25,187) = + reac_rate_local(187) 
  reac_source_local(02,188) = + reac_rate_local(188) 
  reac_source_local(09,188) = - reac_rate_local(188) 
  reac_source_local(17,188) = - reac_rate_local(188) 
  reac_source_local(25,188) = + reac_rate_local(188) 
  reac_source_local(02,189) = + reac_rate_local(189) 
  reac_source_local(09,189) = - reac_rate_local(189) 
  reac_source_local(18,189) = - reac_rate_local(189) 
  reac_source_local(25,189) = + reac_rate_local(189) 
  reac_source_local(02,190) = + reac_rate_local(190) 
  reac_source_local(10,190) = - reac_rate_local(190) 
  reac_source_local(15,190) = - reac_rate_local(190) 
  reac_source_local(25,190) = + reac_rate_local(190) 
  reac_source_local(02,191) = + reac_rate_local(191) 
  reac_source_local(10,191) = - reac_rate_local(191) 
  reac_source_local(16,191) = - reac_rate_local(191) 
  reac_source_local(25,191) = + reac_rate_local(191) 
  reac_source_local(02,192) = + reac_rate_local(192) 
  reac_source_local(10,192) = - reac_rate_local(192) 
  reac_source_local(17,192) = - reac_rate_local(192) 
  reac_source_local(25,192) = + reac_rate_local(192) 
  reac_source_local(02,193) = + reac_rate_local(193) 
  reac_source_local(10,193) = - reac_rate_local(193) 
  reac_source_local(18,193) = - reac_rate_local(193) 
  reac_source_local(25,193) = + reac_rate_local(193) 
  reac_source_local(02,194) = + reac_rate_local(194) 
  reac_source_local(10,194) = - reac_rate_local(194) 
  reac_source_local(15,194) = + reac_rate_local(194) 
  reac_source_local(24,194) = - reac_rate_local(194) 
  reac_source_local(10,195) = - reac_rate_local(195) 
  reac_source_local(21,195) = + reac_rate_local(195) 
  reac_source_local(22,195) = - reac_rate_local(195) 
  reac_source_local(09,196) = + reac_rate_local(196) 
  reac_source_local(15,196) = - reac_rate_local(196) 
  reac_source_local(22,196) = - reac_rate_local(196) 
  reac_source_local(25,196) = + reac_rate_local(196) 
  reac_source_local(02,197) = + reac_rate_local(197) 
  reac_source_local(09,197) = + reac_rate_local(197) 
  reac_source_local(15,197) = - reac_rate_local(197) 
  reac_source_local(21,197) = - reac_rate_local(197) 
  reac_source_local(15,198) = - reac_rate_local(198) 
  reac_source_local(24,198) = + reac_rate_local(198) 
  reac_source_local(25,198) = - reac_rate_local(198) 
  reac_source_local(09,199) = + reac_rate_local(199) 
  reac_source_local(16,199) = - reac_rate_local(199) 
  reac_source_local(22,199) = - reac_rate_local(199) 
  reac_source_local(25,199) = + reac_rate_local(199) 
  reac_source_local(09,200) = + reac_rate_local(200) 
  reac_source_local(17,200) = - reac_rate_local(200) 
  reac_source_local(22,200) = - reac_rate_local(200) 
  reac_source_local(25,200) = + reac_rate_local(200) 
  reac_source_local(09,201) = + reac_rate_local(201) 
  reac_source_local(18,201) = - reac_rate_local(201) 
  reac_source_local(22,201) = - reac_rate_local(201) 
  reac_source_local(25,201) = + reac_rate_local(201) 
  reac_source_local(16,202) = - reac_rate_local(202) 
  reac_source_local(24,202) = + reac_rate_local(202) 
  reac_source_local(25,202) = - reac_rate_local(202) 
  reac_source_local(17,203) = - reac_rate_local(203) 
  reac_source_local(24,203) = + reac_rate_local(203) 
  reac_source_local(25,203) = - reac_rate_local(203) 
  reac_source_local(18,204) = - reac_rate_local(204) 
  reac_source_local(24,204) = + reac_rate_local(204) 
  reac_source_local(25,204) = - reac_rate_local(204) 
  reac_source_local(09,205) = + reac_rate_local(205) 
  reac_source_local(22,205) = - reac_rate_local(205) 
  reac_source_local(25,205) = - reac_rate_local(205) 
  reac_source_local(09,206) = + reac_rate_local(206) * 2.d0
  reac_source_local(21,206) = - reac_rate_local(206) 
  reac_source_local(25,206) = - reac_rate_local(206) 
  reac_source_local(15,207) = + reac_rate_local(207) * 2.d0
  reac_source_local(24,207) = - reac_rate_local(207) 
  reac_source_local(25,207) = - reac_rate_local(207) 
  reac_source_local(15,208) = + reac_rate_local(208) 
  reac_source_local(24,208) = - reac_rate_local(208) 
  reac_source_local(25,208) = + reac_rate_local(208) 
  reac_source_local(15,209) = + reac_rate_local(209) 
  reac_source_local(25,209) = - reac_rate_local(209) * 2.d0
  reac_source_local(02,210) = - reac_rate_local(210) 
  reac_source_local(27,210) = + reac_rate_local(210) 
  reac_source_local(36,210) = - reac_rate_local(210) 
  reac_source_local(02,211) = - reac_rate_local(211) 
  reac_source_local(09,211) = + reac_rate_local(211) 
  reac_source_local(36,211) = + reac_rate_local(211) 
  reac_source_local(37,211) = - reac_rate_local(211) 
  reac_source_local(02,212) = - reac_rate_local(212) 
  reac_source_local(25,212) = + reac_rate_local(212) 
  reac_source_local(26,212) = + reac_rate_local(212) 
  reac_source_local(37,212) = - reac_rate_local(212) 
  reac_source_local(02,213) = - reac_rate_local(213) 
  reac_source_local(09,213) = + reac_rate_local(213) 
  reac_source_local(28,213) = + reac_rate_local(213) 
  reac_source_local(32,213) = - reac_rate_local(213) 
  reac_source_local(02,214) = - reac_rate_local(214) 
  reac_source_local(09,214) = + reac_rate_local(214) 
  reac_source_local(26,214) = + reac_rate_local(214) 
  reac_source_local(28,214) = - reac_rate_local(214) 
  reac_source_local(02,215) = - reac_rate_local(215) 
  reac_source_local(34,215) = + reac_rate_local(215) 
  reac_source_local(39,215) = - reac_rate_local(215) 
  reac_source_local(02,216) = - reac_rate_local(216) 
  reac_source_local(35,216) = + reac_rate_local(216) 
  reac_source_local(40,216) = - reac_rate_local(216) 
  reac_source_local(02,217) = - reac_rate_local(217) 
  reac_source_local(15,217) = + reac_rate_local(217) 
  reac_source_local(34,217) = + reac_rate_local(217) 
  reac_source_local(41,217) = - reac_rate_local(217) 
  reac_source_local(02,218) = - reac_rate_local(218) 
  reac_source_local(15,218) = + reac_rate_local(218) 
  reac_source_local(35,218) = + reac_rate_local(218) 
  reac_source_local(42,218) = - reac_rate_local(218) 
  reac_source_local(02,219) = - reac_rate_local(219) 
  reac_source_local(26,219) = - reac_rate_local(219) 
  reac_source_local(31,219) = + reac_rate_local(219) 
  reac_source_local(03,220) = - reac_rate_local(220) 
  reac_source_local(09,220) = + reac_rate_local(220) 
  reac_source_local(36,220) = + reac_rate_local(220) 
  reac_source_local(37,220) = - reac_rate_local(220) 
  reac_source_local(04,221) = - reac_rate_local(221) 
  reac_source_local(09,221) = + reac_rate_local(221) 
  reac_source_local(36,221) = + reac_rate_local(221) 
  reac_source_local(37,221) = - reac_rate_local(221) 
  reac_source_local(05,222) = - reac_rate_local(222) 
  reac_source_local(09,222) = + reac_rate_local(222) 
  reac_source_local(36,222) = + reac_rate_local(222) 
  reac_source_local(37,222) = - reac_rate_local(222) 
  reac_source_local(06,223) = - reac_rate_local(223) 
  reac_source_local(09,223) = + reac_rate_local(223) 
  reac_source_local(36,223) = + reac_rate_local(223) 
  reac_source_local(37,223) = - reac_rate_local(223) 
  reac_source_local(03,224) = - reac_rate_local(224) 
  reac_source_local(25,224) = + reac_rate_local(224) 
  reac_source_local(26,224) = + reac_rate_local(224) 
  reac_source_local(37,224) = - reac_rate_local(224) 
  reac_source_local(04,225) = - reac_rate_local(225) 
  reac_source_local(25,225) = + reac_rate_local(225) 
  reac_source_local(26,225) = + reac_rate_local(225) 
  reac_source_local(37,225) = - reac_rate_local(225) 
  reac_source_local(05,226) = - reac_rate_local(226) 
  reac_source_local(25,226) = + reac_rate_local(226) 
  reac_source_local(26,226) = + reac_rate_local(226) 
  reac_source_local(37,226) = - reac_rate_local(226) 
  reac_source_local(06,227) = - reac_rate_local(227) 
  reac_source_local(25,227) = + reac_rate_local(227) 
  reac_source_local(26,227) = + reac_rate_local(227) 
  reac_source_local(37,227) = - reac_rate_local(227) 
  reac_source_local(03,228) = - reac_rate_local(228) 
  reac_source_local(09,228) = + reac_rate_local(228) 
  reac_source_local(28,228) = + reac_rate_local(228) 
  reac_source_local(32,228) = - reac_rate_local(228) 
  reac_source_local(04,229) = - reac_rate_local(229) 
  reac_source_local(09,229) = + reac_rate_local(229) 
  reac_source_local(28,229) = + reac_rate_local(229) 
  reac_source_local(32,229) = - reac_rate_local(229) 
  reac_source_local(05,230) = - reac_rate_local(230) 
  reac_source_local(09,230) = + reac_rate_local(230) 
  reac_source_local(28,230) = + reac_rate_local(230) 
  reac_source_local(32,230) = - reac_rate_local(230) 
  reac_source_local(06,231) = - reac_rate_local(231) 
  reac_source_local(09,231) = + reac_rate_local(231) 
  reac_source_local(28,231) = + reac_rate_local(231) 
  reac_source_local(32,231) = - reac_rate_local(231) 
  reac_source_local(03,232) = - reac_rate_local(232) 
  reac_source_local(09,232) = + reac_rate_local(232) 
  reac_source_local(26,232) = + reac_rate_local(232) 
  reac_source_local(28,232) = - reac_rate_local(232) 
  reac_source_local(04,233) = - reac_rate_local(233) 
  reac_source_local(09,233) = + reac_rate_local(233) 
  reac_source_local(26,233) = + reac_rate_local(233) 
  reac_source_local(28,233) = - reac_rate_local(233) 
  reac_source_local(05,234) = - reac_rate_local(234) 
  reac_source_local(09,234) = + reac_rate_local(234) 
  reac_source_local(26,234) = + reac_rate_local(234) 
  reac_source_local(28,234) = - reac_rate_local(234) 
  reac_source_local(06,235) = - reac_rate_local(235) 
  reac_source_local(09,235) = + reac_rate_local(235) 
  reac_source_local(26,235) = + reac_rate_local(235) 
  reac_source_local(28,235) = - reac_rate_local(235) 
  reac_source_local(03,236) = - reac_rate_local(236) 
  reac_source_local(34,236) = + reac_rate_local(236) 
  reac_source_local(39,236) = - reac_rate_local(236) 
  reac_source_local(04,237) = - reac_rate_local(237) 
  reac_source_local(34,237) = + reac_rate_local(237) 
  reac_source_local(39,237) = - reac_rate_local(237) 
  reac_source_local(05,238) = - reac_rate_local(238) 
  reac_source_local(34,238) = + reac_rate_local(238) 
  reac_source_local(39,238) = - reac_rate_local(238) 
  reac_source_local(06,239) = - reac_rate_local(239) 
  reac_source_local(34,239) = + reac_rate_local(239) 
  reac_source_local(39,239) = - reac_rate_local(239) 
  reac_source_local(03,240) = - reac_rate_local(240) 
  reac_source_local(35,240) = + reac_rate_local(240) 
  reac_source_local(40,240) = - reac_rate_local(240) 
  reac_source_local(04,241) = - reac_rate_local(241) 
  reac_source_local(35,241) = + reac_rate_local(241) 
  reac_source_local(40,241) = - reac_rate_local(241) 
  reac_source_local(05,242) = - reac_rate_local(242) 
  reac_source_local(35,242) = + reac_rate_local(242) 
  reac_source_local(40,242) = - reac_rate_local(242) 
  reac_source_local(06,243) = - reac_rate_local(243) 
  reac_source_local(35,243) = + reac_rate_local(243) 
  reac_source_local(40,243) = - reac_rate_local(243) 
  reac_source_local(03,244) = - reac_rate_local(244) 
  reac_source_local(15,244) = + reac_rate_local(244) 
  reac_source_local(34,244) = + reac_rate_local(244) 
  reac_source_local(41,244) = - reac_rate_local(244) 
  reac_source_local(04,245) = - reac_rate_local(245) 
  reac_source_local(15,245) = + reac_rate_local(245) 
  reac_source_local(34,245) = + reac_rate_local(245) 
  reac_source_local(41,245) = - reac_rate_local(245) 
  reac_source_local(05,246) = - reac_rate_local(246) 
  reac_source_local(15,246) = + reac_rate_local(246) 
  reac_source_local(34,246) = + reac_rate_local(246) 
  reac_source_local(41,246) = - reac_rate_local(246) 
  reac_source_local(06,247) = - reac_rate_local(247) 
  reac_source_local(15,247) = + reac_rate_local(247) 
  reac_source_local(34,247) = + reac_rate_local(247) 
  reac_source_local(41,247) = - reac_rate_local(247) 
  reac_source_local(03,248) = - reac_rate_local(248) 
  reac_source_local(15,248) = + reac_rate_local(248) 
  reac_source_local(35,248) = + reac_rate_local(248) 
  reac_source_local(42,248) = - reac_rate_local(248) 
  reac_source_local(04,249) = - reac_rate_local(249) 
  reac_source_local(15,249) = + reac_rate_local(249) 
  reac_source_local(35,249) = + reac_rate_local(249) 
  reac_source_local(42,249) = - reac_rate_local(249) 
  reac_source_local(05,250) = - reac_rate_local(250) 
  reac_source_local(15,250) = + reac_rate_local(250) 
  reac_source_local(35,250) = + reac_rate_local(250) 
  reac_source_local(42,250) = - reac_rate_local(250) 
  reac_source_local(06,251) = - reac_rate_local(251) 
  reac_source_local(15,251) = + reac_rate_local(251) 
  reac_source_local(35,251) = + reac_rate_local(251) 
  reac_source_local(42,251) = - reac_rate_local(251) 
  reac_source_local(03,252) = - reac_rate_local(252) 
  reac_source_local(26,252) = - reac_rate_local(252) 
  reac_source_local(31,252) = + reac_rate_local(252) 
  reac_source_local(04,253) = - reac_rate_local(253) 
  reac_source_local(26,253) = - reac_rate_local(253) 
  reac_source_local(31,253) = + reac_rate_local(253) 
  reac_source_local(05,254) = - reac_rate_local(254) 
  reac_source_local(26,254) = - reac_rate_local(254) 
  reac_source_local(31,254) = + reac_rate_local(254) 
  reac_source_local(06,255) = - reac_rate_local(255) 
  reac_source_local(26,255) = - reac_rate_local(255) 
  reac_source_local(31,255) = + reac_rate_local(255) 
  reac_source_local(03,256) = - reac_rate_local(256) 
  reac_source_local(27,256) = + reac_rate_local(256) 
  reac_source_local(36,256) = - reac_rate_local(256) 
  reac_source_local(04,257) = - reac_rate_local(257) 
  reac_source_local(27,257) = + reac_rate_local(257) 
  reac_source_local(36,257) = - reac_rate_local(257) 
  reac_source_local(05,258) = - reac_rate_local(258) 
  reac_source_local(27,258) = + reac_rate_local(258) 
  reac_source_local(36,258) = - reac_rate_local(258) 
  reac_source_local(06,259) = - reac_rate_local(259) 
  reac_source_local(27,259) = + reac_rate_local(259) 
  reac_source_local(36,259) = - reac_rate_local(259) 
  reac_source_local(07,260) = - reac_rate_local(260) 
  reac_source_local(09,260) = + reac_rate_local(260) 
  reac_source_local(36,260) = + reac_rate_local(260) 
  reac_source_local(37,260) = - reac_rate_local(260) 
  reac_source_local(08,261) = - reac_rate_local(261) 
  reac_source_local(09,261) = + reac_rate_local(261) 
  reac_source_local(36,261) = + reac_rate_local(261) 
  reac_source_local(37,261) = - reac_rate_local(261) 
  reac_source_local(07,262) = - reac_rate_local(262) 
  reac_source_local(25,262) = + reac_rate_local(262) 
  reac_source_local(26,262) = + reac_rate_local(262) 
  reac_source_local(37,262) = - reac_rate_local(262) 
  reac_source_local(08,263) = - reac_rate_local(263) 
  reac_source_local(25,263) = + reac_rate_local(263) 
  reac_source_local(26,263) = + reac_rate_local(263) 
  reac_source_local(37,263) = - reac_rate_local(263) 
  reac_source_local(07,264) = - reac_rate_local(264) 
  reac_source_local(09,264) = + reac_rate_local(264) 
  reac_source_local(28,264) = + reac_rate_local(264) 
  reac_source_local(32,264) = - reac_rate_local(264) 
  reac_source_local(08,265) = - reac_rate_local(265) 
  reac_source_local(09,265) = + reac_rate_local(265) 
  reac_source_local(28,265) = + reac_rate_local(265) 
  reac_source_local(32,265) = - reac_rate_local(265) 
  reac_source_local(07,266) = - reac_rate_local(266) 
  reac_source_local(09,266) = + reac_rate_local(266) 
  reac_source_local(26,266) = + reac_rate_local(266) 
  reac_source_local(28,266) = - reac_rate_local(266) 
  reac_source_local(08,267) = - reac_rate_local(267) 
  reac_source_local(09,267) = + reac_rate_local(267) 
  reac_source_local(26,267) = + reac_rate_local(267) 
  reac_source_local(28,267) = - reac_rate_local(267) 
  reac_source_local(07,268) = - reac_rate_local(268) 
  reac_source_local(34,268) = + reac_rate_local(268) 
  reac_source_local(39,268) = - reac_rate_local(268) 
  reac_source_local(08,269) = - reac_rate_local(269) 
  reac_source_local(34,269) = + reac_rate_local(269) 
  reac_source_local(39,269) = - reac_rate_local(269) 
  reac_source_local(07,270) = - reac_rate_local(270) 
  reac_source_local(35,270) = + reac_rate_local(270) 
  reac_source_local(40,270) = - reac_rate_local(270) 
  reac_source_local(08,271) = - reac_rate_local(271) 
  reac_source_local(35,271) = + reac_rate_local(271) 
  reac_source_local(40,271) = - reac_rate_local(271) 
  reac_source_local(07,272) = - reac_rate_local(272) 
  reac_source_local(15,272) = + reac_rate_local(272) 
  reac_source_local(34,272) = + reac_rate_local(272) 
  reac_source_local(41,272) = - reac_rate_local(272) 
  reac_source_local(08,273) = - reac_rate_local(273) 
  reac_source_local(15,273) = + reac_rate_local(273) 
  reac_source_local(34,273) = + reac_rate_local(273) 
  reac_source_local(41,273) = - reac_rate_local(273) 
  reac_source_local(07,274) = - reac_rate_local(274) 
  reac_source_local(15,274) = + reac_rate_local(274) 
  reac_source_local(35,274) = + reac_rate_local(274) 
  reac_source_local(42,274) = - reac_rate_local(274) 
  reac_source_local(08,275) = - reac_rate_local(275) 
  reac_source_local(15,275) = + reac_rate_local(275) 
  reac_source_local(35,275) = + reac_rate_local(275) 
  reac_source_local(42,275) = - reac_rate_local(275) 
  reac_source_local(07,276) = - reac_rate_local(276) 
  reac_source_local(26,276) = - reac_rate_local(276) 
  reac_source_local(31,276) = + reac_rate_local(276) 
  reac_source_local(08,277) = - reac_rate_local(277) 
  reac_source_local(26,277) = - reac_rate_local(277) 
  reac_source_local(31,277) = + reac_rate_local(277) 
  reac_source_local(07,278) = - reac_rate_local(278) 
  reac_source_local(27,278) = + reac_rate_local(278) 
  reac_source_local(36,278) = - reac_rate_local(278) 
  reac_source_local(08,279) = - reac_rate_local(279) 
  reac_source_local(27,279) = + reac_rate_local(279) 
  reac_source_local(36,279) = - reac_rate_local(279) 
  reac_source_local(07,280) = - reac_rate_local(280) 
  reac_source_local(27,280) = + reac_rate_local(280) 
  reac_source_local(36,280) = - reac_rate_local(280) 
  reac_source_local(08,281) = - reac_rate_local(281) 
  reac_source_local(27,281) = + reac_rate_local(281) 
  reac_source_local(36,281) = - reac_rate_local(281) 
  reac_source_local(09,282) = - reac_rate_local(282) 
  reac_source_local(25,282) = + reac_rate_local(282) 
  reac_source_local(28,282) = + reac_rate_local(282) 
  reac_source_local(37,282) = - reac_rate_local(282) 
  reac_source_local(01,283) = + reac_rate_local(283) 
  reac_source_local(02,283) = + reac_rate_local(283) 
  reac_source_local(09,283) = - reac_rate_local(283) 
  reac_source_local(39,283) = - reac_rate_local(283) 
  reac_source_local(01,284) = + reac_rate_local(284) 
  reac_source_local(02,284) = + reac_rate_local(284) * 2.d0
  reac_source_local(09,284) = - reac_rate_local(284) 
  reac_source_local(34,284) = - reac_rate_local(284) 
  reac_source_local(02,285) = + reac_rate_local(285) 
  reac_source_local(09,285) = - reac_rate_local(285) 
  reac_source_local(29,285) = + reac_rate_local(285) 
  reac_source_local(30,285) = - reac_rate_local(285) 
  reac_source_local(02,286) = + reac_rate_local(286) 
  reac_source_local(09,286) = - reac_rate_local(286) 
  reac_source_local(30,286) = + reac_rate_local(286) 
  reac_source_local(31,286) = - reac_rate_local(286) 
  reac_source_local(02,287) = + reac_rate_local(287) 
  reac_source_local(09,287) = - reac_rate_local(287) 
  reac_source_local(29,287) = + reac_rate_local(287) 
  reac_source_local(30,287) = - reac_rate_local(287) 
  reac_source_local(02,288) = + reac_rate_local(288) 
  reac_source_local(09,288) = - reac_rate_local(288) 
  reac_source_local(30,288) = + reac_rate_local(288) 
  reac_source_local(31,288) = - reac_rate_local(288) 
  reac_source_local(09,289) = - reac_rate_local(289) 
  reac_source_local(22,289) = + reac_rate_local(289) 
  reac_source_local(28,289) = + reac_rate_local(289) 
  reac_source_local(32,289) = - reac_rate_local(289) 
  reac_source_local(10,290) = - reac_rate_local(290) 
  reac_source_local(25,290) = + reac_rate_local(290) 
  reac_source_local(28,290) = + reac_rate_local(290) 
  reac_source_local(37,290) = - reac_rate_local(290) 
  reac_source_local(01,291) = + reac_rate_local(291) 
  reac_source_local(02,291) = + reac_rate_local(291) 
  reac_source_local(10,291) = - reac_rate_local(291) 
  reac_source_local(39,291) = - reac_rate_local(291) 
  reac_source_local(01,292) = + reac_rate_local(292) 
  reac_source_local(02,292) = + reac_rate_local(292) * 2.d0
  reac_source_local(10,292) = - reac_rate_local(292) 
  reac_source_local(34,292) = - reac_rate_local(292) 
  reac_source_local(02,293) = + reac_rate_local(293) 
  reac_source_local(10,293) = - reac_rate_local(293) 
  reac_source_local(29,293) = + reac_rate_local(293) 
  reac_source_local(30,293) = - reac_rate_local(293) 
  reac_source_local(02,294) = + reac_rate_local(294) 
  reac_source_local(10,294) = - reac_rate_local(294) 
  reac_source_local(30,294) = + reac_rate_local(294) 
  reac_source_local(31,294) = - reac_rate_local(294) 
  reac_source_local(02,295) = + reac_rate_local(295) 
  reac_source_local(10,295) = - reac_rate_local(295) 
  reac_source_local(29,295) = + reac_rate_local(295) 
  reac_source_local(30,295) = - reac_rate_local(295) 
  reac_source_local(02,296) = + reac_rate_local(296) 
  reac_source_local(10,296) = - reac_rate_local(296) 
  reac_source_local(30,296) = + reac_rate_local(296) 
  reac_source_local(31,296) = - reac_rate_local(296) 
  reac_source_local(10,297) = - reac_rate_local(297) 
  reac_source_local(22,297) = + reac_rate_local(297) 
  reac_source_local(28,297) = + reac_rate_local(297) 
  reac_source_local(32,297) = - reac_rate_local(297) 
  reac_source_local(01,298) = + reac_rate_local(298) 
  reac_source_local(02,298) = + reac_rate_local(298) 
  reac_source_local(11,298) = - reac_rate_local(298) 
  reac_source_local(39,298) = - reac_rate_local(298) 
  reac_source_local(01,299) = + reac_rate_local(299) 
  reac_source_local(02,299) = + reac_rate_local(299) 
  reac_source_local(12,299) = - reac_rate_local(299) 
  reac_source_local(39,299) = - reac_rate_local(299) 
  reac_source_local(01,300) = + reac_rate_local(300) 
  reac_source_local(02,300) = + reac_rate_local(300) 
  reac_source_local(13,300) = - reac_rate_local(300) 
  reac_source_local(39,300) = - reac_rate_local(300) 
  reac_source_local(01,301) = + reac_rate_local(301) 
  reac_source_local(02,301) = + reac_rate_local(301) 
  reac_source_local(14,301) = - reac_rate_local(301) 
  reac_source_local(39,301) = - reac_rate_local(301) 
  reac_source_local(01,302) = + reac_rate_local(302) 
  reac_source_local(02,302) = + reac_rate_local(302) * 2.d0
  reac_source_local(11,302) = - reac_rate_local(302) 
  reac_source_local(34,302) = - reac_rate_local(302) 
  reac_source_local(01,303) = + reac_rate_local(303) 
  reac_source_local(02,303) = + reac_rate_local(303) * 2.d0
  reac_source_local(12,303) = - reac_rate_local(303) 
  reac_source_local(34,303) = - reac_rate_local(303) 
  reac_source_local(01,304) = + reac_rate_local(304) 
  reac_source_local(02,304) = + reac_rate_local(304) * 2.d0
  reac_source_local(13,304) = - reac_rate_local(304) 
  reac_source_local(34,304) = - reac_rate_local(304) 
  reac_source_local(01,305) = + reac_rate_local(305) 
  reac_source_local(02,305) = + reac_rate_local(305) * 2.d0
  reac_source_local(14,305) = - reac_rate_local(305) 
  reac_source_local(34,305) = - reac_rate_local(305) 
  reac_source_local(02,306) = + reac_rate_local(306) 
  reac_source_local(11,306) = - reac_rate_local(306) 
  reac_source_local(29,306) = + reac_rate_local(306) 
  reac_source_local(30,306) = - reac_rate_local(306) 
  reac_source_local(02,307) = + reac_rate_local(307) 
  reac_source_local(12,307) = - reac_rate_local(307) 
  reac_source_local(29,307) = + reac_rate_local(307) 
  reac_source_local(30,307) = - reac_rate_local(307) 
  reac_source_local(02,308) = + reac_rate_local(308) 
  reac_source_local(13,308) = - reac_rate_local(308) 
  reac_source_local(29,308) = + reac_rate_local(308) 
  reac_source_local(30,308) = - reac_rate_local(308) 
  reac_source_local(02,309) = + reac_rate_local(309) 
  reac_source_local(14,309) = - reac_rate_local(309) 
  reac_source_local(29,309) = + reac_rate_local(309) 
  reac_source_local(30,309) = - reac_rate_local(309) 
  reac_source_local(02,310) = + reac_rate_local(310) 
  reac_source_local(11,310) = - reac_rate_local(310) 
  reac_source_local(30,310) = + reac_rate_local(310) 
  reac_source_local(31,310) = - reac_rate_local(310) 
  reac_source_local(02,311) = + reac_rate_local(311) 
  reac_source_local(12,311) = - reac_rate_local(311) 
  reac_source_local(30,311) = + reac_rate_local(311) 
  reac_source_local(31,311) = - reac_rate_local(311) 
  reac_source_local(02,312) = + reac_rate_local(312) 
  reac_source_local(13,312) = - reac_rate_local(312) 
  reac_source_local(30,312) = + reac_rate_local(312) 
  reac_source_local(31,312) = - reac_rate_local(312) 
  reac_source_local(02,313) = + reac_rate_local(313) 
  reac_source_local(14,313) = - reac_rate_local(313) 
  reac_source_local(30,313) = + reac_rate_local(313) 
  reac_source_local(31,313) = - reac_rate_local(313) 
  reac_source_local(02,314) = + reac_rate_local(314) 
  reac_source_local(11,314) = - reac_rate_local(314) 
  reac_source_local(29,314) = + reac_rate_local(314) 
  reac_source_local(30,314) = - reac_rate_local(314) 
  reac_source_local(02,315) = + reac_rate_local(315) 
  reac_source_local(12,315) = - reac_rate_local(315) 
  reac_source_local(29,315) = + reac_rate_local(315) 
  reac_source_local(30,315) = - reac_rate_local(315) 
  reac_source_local(02,316) = + reac_rate_local(316) 
  reac_source_local(13,316) = - reac_rate_local(316) 
  reac_source_local(29,316) = + reac_rate_local(316) 
  reac_source_local(30,316) = - reac_rate_local(316) 
  reac_source_local(02,317) = + reac_rate_local(317) 
  reac_source_local(14,317) = - reac_rate_local(317) 
  reac_source_local(29,317) = + reac_rate_local(317) 
  reac_source_local(30,317) = - reac_rate_local(317) 
  reac_source_local(02,318) = + reac_rate_local(318) 
  reac_source_local(11,318) = - reac_rate_local(318) 
  reac_source_local(30,318) = + reac_rate_local(318) 
  reac_source_local(31,318) = - reac_rate_local(318) 
  reac_source_local(02,319) = + reac_rate_local(319) 
  reac_source_local(12,319) = - reac_rate_local(319) 
  reac_source_local(30,319) = + reac_rate_local(319) 
  reac_source_local(31,319) = - reac_rate_local(319) 
  reac_source_local(02,320) = + reac_rate_local(320) 
  reac_source_local(13,320) = - reac_rate_local(320) 
  reac_source_local(30,320) = + reac_rate_local(320) 
  reac_source_local(31,320) = - reac_rate_local(320) 
  reac_source_local(02,321) = + reac_rate_local(321) 
  reac_source_local(14,321) = - reac_rate_local(321) 
  reac_source_local(30,321) = + reac_rate_local(321) 
  reac_source_local(31,321) = - reac_rate_local(321) 
  reac_source_local(11,322) = - reac_rate_local(322) 
  reac_source_local(22,322) = + reac_rate_local(322) 
  reac_source_local(28,322) = + reac_rate_local(322) 
  reac_source_local(32,322) = - reac_rate_local(322) 
  reac_source_local(12,323) = - reac_rate_local(323) 
  reac_source_local(22,323) = + reac_rate_local(323) 
  reac_source_local(28,323) = + reac_rate_local(323) 
  reac_source_local(32,323) = - reac_rate_local(323) 
  reac_source_local(13,324) = - reac_rate_local(324) 
  reac_source_local(22,324) = + reac_rate_local(324) 
  reac_source_local(28,324) = + reac_rate_local(324) 
  reac_source_local(32,324) = - reac_rate_local(324) 
  reac_source_local(14,325) = - reac_rate_local(325) 
  reac_source_local(22,325) = + reac_rate_local(325) 
  reac_source_local(28,325) = + reac_rate_local(325) 
  reac_source_local(32,325) = - reac_rate_local(325) 
  reac_source_local(02,326) = + reac_rate_local(326) 
  reac_source_local(15,326) = - reac_rate_local(326) 
  reac_source_local(26,326) = - reac_rate_local(326) 
  reac_source_local(36,326) = + reac_rate_local(326) 
  reac_source_local(09,327) = + reac_rate_local(327) 
  reac_source_local(15,327) = - reac_rate_local(327) 
  reac_source_local(28,327) = - reac_rate_local(327) 
  reac_source_local(36,327) = + reac_rate_local(327) 
  reac_source_local(09,328) = + reac_rate_local(328) * 2.d0
  reac_source_local(15,328) = - reac_rate_local(328) 
  reac_source_local(29,328) = - reac_rate_local(328) 
  reac_source_local(36,328) = + reac_rate_local(328) 
  reac_source_local(09,329) = + reac_rate_local(329) 
  reac_source_local(15,329) = - reac_rate_local(329) 
  reac_source_local(32,329) = - reac_rate_local(329) 
  reac_source_local(37,329) = + reac_rate_local(329) 
  reac_source_local(15,330) = - reac_rate_local(330) 
  reac_source_local(25,330) = + reac_rate_local(330) 
  reac_source_local(28,330) = + reac_rate_local(330) 
  reac_source_local(32,330) = - reac_rate_local(330) 
  reac_source_local(15,331) = - reac_rate_local(331) 
  reac_source_local(25,331) = + reac_rate_local(331) 
  reac_source_local(36,331) = + reac_rate_local(331) 
  reac_source_local(37,331) = - reac_rate_local(331) 
  reac_source_local(15,332) = - reac_rate_local(332) 
  reac_source_local(36,332) = - reac_rate_local(332) 
  reac_source_local(38,332) = + reac_rate_local(332) 
  reac_source_local(15,333) = - reac_rate_local(333) 
  reac_source_local(40,333) = - reac_rate_local(333) 
  reac_source_local(42,333) = + reac_rate_local(333) 
  reac_source_local(01,334) = + reac_rate_local(334) 
  reac_source_local(15,334) = - reac_rate_local(334) 
  reac_source_local(24,334) = + reac_rate_local(334) 
  reac_source_local(39,334) = - reac_rate_local(334) 
  reac_source_local(15,335) = - reac_rate_local(335) 
  reac_source_local(39,335) = - reac_rate_local(335) 
  reac_source_local(41,335) = + reac_rate_local(335) 
  reac_source_local(01,336) = + reac_rate_local(336) 
  reac_source_local(15,336) = + reac_rate_local(336) 
  reac_source_local(40,336) = - reac_rate_local(336) 
  reac_source_local(01,337) = + reac_rate_local(337) 
  reac_source_local(24,337) = + reac_rate_local(337) 
  reac_source_local(41,337) = - reac_rate_local(337) 
  reac_source_local(02,338) = + reac_rate_local(338) 
  reac_source_local(16,338) = - reac_rate_local(338) 
  reac_source_local(26,338) = - reac_rate_local(338) 
  reac_source_local(36,338) = + reac_rate_local(338) 
  reac_source_local(02,339) = + reac_rate_local(339) 
  reac_source_local(17,339) = - reac_rate_local(339) 
  reac_source_local(26,339) = - reac_rate_local(339) 
  reac_source_local(36,339) = + reac_rate_local(339) 
  reac_source_local(02,340) = + reac_rate_local(340) 
  reac_source_local(18,340) = - reac_rate_local(340) 
  reac_source_local(26,340) = - reac_rate_local(340) 
  reac_source_local(36,340) = + reac_rate_local(340) 
  reac_source_local(09,341) = + reac_rate_local(341) 
  reac_source_local(16,341) = - reac_rate_local(341) 
  reac_source_local(28,341) = - reac_rate_local(341) 
  reac_source_local(36,341) = + reac_rate_local(341) 
  reac_source_local(09,342) = + reac_rate_local(342) 
  reac_source_local(17,342) = - reac_rate_local(342) 
  reac_source_local(28,342) = - reac_rate_local(342) 
  reac_source_local(36,342) = + reac_rate_local(342) 
  reac_source_local(09,343) = + reac_rate_local(343) 
  reac_source_local(18,343) = - reac_rate_local(343) 
  reac_source_local(28,343) = - reac_rate_local(343) 
  reac_source_local(36,343) = + reac_rate_local(343) 
  reac_source_local(09,344) = + reac_rate_local(344) * 2.d0
  reac_source_local(16,344) = - reac_rate_local(344) 
  reac_source_local(29,344) = - reac_rate_local(344) 
  reac_source_local(36,344) = + reac_rate_local(344) 
  reac_source_local(09,345) = + reac_rate_local(345) * 2.d0
  reac_source_local(17,345) = - reac_rate_local(345) 
  reac_source_local(29,345) = - reac_rate_local(345) 
  reac_source_local(36,345) = + reac_rate_local(345) 
  reac_source_local(09,346) = + reac_rate_local(346) * 2.d0
  reac_source_local(18,346) = - reac_rate_local(346) 
  reac_source_local(29,346) = - reac_rate_local(346) 
  reac_source_local(36,346) = + reac_rate_local(346) 
  reac_source_local(09,347) = + reac_rate_local(347) 
  reac_source_local(16,347) = - reac_rate_local(347) 
  reac_source_local(32,347) = - reac_rate_local(347) 
  reac_source_local(37,347) = + reac_rate_local(347) 
  reac_source_local(09,348) = + reac_rate_local(348) 
  reac_source_local(17,348) = - reac_rate_local(348) 
  reac_source_local(32,348) = - reac_rate_local(348) 
  reac_source_local(37,348) = + reac_rate_local(348) 
  reac_source_local(09,349) = + reac_rate_local(349) 
  reac_source_local(18,349) = - reac_rate_local(349) 
  reac_source_local(32,349) = - reac_rate_local(349) 
  reac_source_local(37,349) = + reac_rate_local(349) 
  reac_source_local(16,350) = - reac_rate_local(350) 
  reac_source_local(25,350) = + reac_rate_local(350) 
  reac_source_local(28,350) = + reac_rate_local(350) 
  reac_source_local(32,350) = - reac_rate_local(350) 
  reac_source_local(17,351) = - reac_rate_local(351) 
  reac_source_local(25,351) = + reac_rate_local(351) 
  reac_source_local(28,351) = + reac_rate_local(351) 
  reac_source_local(32,351) = - reac_rate_local(351) 
  reac_source_local(18,352) = - reac_rate_local(352) 
  reac_source_local(25,352) = + reac_rate_local(352) 
  reac_source_local(28,352) = + reac_rate_local(352) 
  reac_source_local(32,352) = - reac_rate_local(352) 
  reac_source_local(16,353) = - reac_rate_local(353) 
  reac_source_local(25,353) = + reac_rate_local(353) 
  reac_source_local(36,353) = + reac_rate_local(353) 
  reac_source_local(37,353) = - reac_rate_local(353) 
  reac_source_local(17,354) = - reac_rate_local(354) 
  reac_source_local(25,354) = + reac_rate_local(354) 
  reac_source_local(36,354) = + reac_rate_local(354) 
  reac_source_local(37,354) = - reac_rate_local(354) 
  reac_source_local(18,355) = - reac_rate_local(355) 
  reac_source_local(25,355) = + reac_rate_local(355) 
  reac_source_local(36,355) = + reac_rate_local(355) 
  reac_source_local(37,355) = - reac_rate_local(355) 
  reac_source_local(16,356) = - reac_rate_local(356) 
  reac_source_local(36,356) = - reac_rate_local(356) 
  reac_source_local(38,356) = + reac_rate_local(356) 
  reac_source_local(17,357) = - reac_rate_local(357) 
  reac_source_local(36,357) = - reac_rate_local(357) 
  reac_source_local(38,357) = + reac_rate_local(357) 
  reac_source_local(18,358) = - reac_rate_local(358) 
  reac_source_local(36,358) = - reac_rate_local(358) 
  reac_source_local(38,358) = + reac_rate_local(358) 
  reac_source_local(16,359) = - reac_rate_local(359) 
  reac_source_local(40,359) = - reac_rate_local(359) 
  reac_source_local(42,359) = + reac_rate_local(359) 
  reac_source_local(17,360) = - reac_rate_local(360) 
  reac_source_local(40,360) = - reac_rate_local(360) 
  reac_source_local(42,360) = + reac_rate_local(360) 
  reac_source_local(18,361) = - reac_rate_local(361) 
  reac_source_local(40,361) = - reac_rate_local(361) 
  reac_source_local(42,361) = + reac_rate_local(361) 
  reac_source_local(01,362) = + reac_rate_local(362) 
  reac_source_local(16,362) = - reac_rate_local(362) 
  reac_source_local(24,362) = + reac_rate_local(362) 
  reac_source_local(39,362) = - reac_rate_local(362) 
  reac_source_local(01,363) = + reac_rate_local(363) 
  reac_source_local(17,363) = - reac_rate_local(363) 
  reac_source_local(24,363) = + reac_rate_local(363) 
  reac_source_local(39,363) = - reac_rate_local(363) 
  reac_source_local(01,364) = + reac_rate_local(364) 
  reac_source_local(18,364) = - reac_rate_local(364) 
  reac_source_local(24,364) = + reac_rate_local(364) 
  reac_source_local(39,364) = - reac_rate_local(364) 
  reac_source_local(16,365) = - reac_rate_local(365) 
  reac_source_local(39,365) = - reac_rate_local(365) 
  reac_source_local(41,365) = + reac_rate_local(365) 
  reac_source_local(17,366) = - reac_rate_local(366) 
  reac_source_local(39,366) = - reac_rate_local(366) 
  reac_source_local(41,366) = + reac_rate_local(366) 
  reac_source_local(18,367) = - reac_rate_local(367) 
  reac_source_local(39,367) = - reac_rate_local(367) 
  reac_source_local(41,367) = + reac_rate_local(367) 
  reac_source_local(01,368) = + reac_rate_local(368) 
  reac_source_local(15,368) = + reac_rate_local(368) 
  reac_source_local(40,368) = - reac_rate_local(368) 
  reac_source_local(01,369) = + reac_rate_local(369) 
  reac_source_local(15,369) = + reac_rate_local(369) 
  reac_source_local(40,369) = - reac_rate_local(369) 
  reac_source_local(01,370) = + reac_rate_local(370) 
  reac_source_local(15,370) = + reac_rate_local(370) 
  reac_source_local(40,370) = - reac_rate_local(370) 
  reac_source_local(01,371) = + reac_rate_local(371) 
  reac_source_local(24,371) = + reac_rate_local(371) 
  reac_source_local(41,371) = - reac_rate_local(371) 
  reac_source_local(01,372) = + reac_rate_local(372) 
  reac_source_local(24,372) = + reac_rate_local(372) 
  reac_source_local(41,372) = - reac_rate_local(372) 
  reac_source_local(01,373) = + reac_rate_local(373) 
  reac_source_local(24,373) = + reac_rate_local(373) 
  reac_source_local(41,373) = - reac_rate_local(373) 
  reac_source_local(02,374) = + reac_rate_local(374) 
  reac_source_local(19,374) = - reac_rate_local(374) 
  reac_source_local(26,374) = - reac_rate_local(374) 
  reac_source_local(36,374) = + reac_rate_local(374) 
  reac_source_local(02,375) = + reac_rate_local(375) 
  reac_source_local(20,375) = - reac_rate_local(375) 
  reac_source_local(26,375) = - reac_rate_local(375) 
  reac_source_local(36,375) = + reac_rate_local(375) 
  reac_source_local(09,376) = + reac_rate_local(376) 
  reac_source_local(19,376) = - reac_rate_local(376) 
  reac_source_local(28,376) = - reac_rate_local(376) 
  reac_source_local(36,376) = + reac_rate_local(376) 
  reac_source_local(09,377) = + reac_rate_local(377) 
  reac_source_local(20,377) = - reac_rate_local(377) 
  reac_source_local(28,377) = - reac_rate_local(377) 
  reac_source_local(36,377) = + reac_rate_local(377) 
  reac_source_local(09,378) = + reac_rate_local(378) * 2.d0
  reac_source_local(19,378) = - reac_rate_local(378) 
  reac_source_local(29,378) = - reac_rate_local(378) 
  reac_source_local(36,378) = + reac_rate_local(378) 
  reac_source_local(09,379) = + reac_rate_local(379) * 2.d0
  reac_source_local(20,379) = - reac_rate_local(379) 
  reac_source_local(29,379) = - reac_rate_local(379) 
  reac_source_local(36,379) = + reac_rate_local(379) 
  reac_source_local(09,380) = + reac_rate_local(380) 
  reac_source_local(19,380) = - reac_rate_local(380) 
  reac_source_local(32,380) = - reac_rate_local(380) 
  reac_source_local(37,380) = + reac_rate_local(380) 
  reac_source_local(09,381) = + reac_rate_local(381) 
  reac_source_local(20,381) = - reac_rate_local(381) 
  reac_source_local(32,381) = - reac_rate_local(381) 
  reac_source_local(37,381) = + reac_rate_local(381) 
  reac_source_local(19,382) = - reac_rate_local(382) 
  reac_source_local(25,382) = + reac_rate_local(382) 
  reac_source_local(28,382) = + reac_rate_local(382) 
  reac_source_local(32,382) = - reac_rate_local(382) 
  reac_source_local(20,383) = - reac_rate_local(383) 
  reac_source_local(25,383) = + reac_rate_local(383) 
  reac_source_local(28,383) = + reac_rate_local(383) 
  reac_source_local(32,383) = - reac_rate_local(383) 
  reac_source_local(19,384) = - reac_rate_local(384) 
  reac_source_local(25,384) = + reac_rate_local(384) 
  reac_source_local(36,384) = + reac_rate_local(384) 
  reac_source_local(37,384) = - reac_rate_local(384) 
  reac_source_local(20,385) = - reac_rate_local(385) 
  reac_source_local(25,385) = + reac_rate_local(385) 
  reac_source_local(36,385) = + reac_rate_local(385) 
  reac_source_local(37,385) = - reac_rate_local(385) 
  reac_source_local(19,386) = - reac_rate_local(386) 
  reac_source_local(36,386) = - reac_rate_local(386) 
  reac_source_local(38,386) = + reac_rate_local(386) 
  reac_source_local(20,387) = - reac_rate_local(387) 
  reac_source_local(36,387) = - reac_rate_local(387) 
  reac_source_local(38,387) = + reac_rate_local(387) 
  reac_source_local(19,388) = - reac_rate_local(388) 
  reac_source_local(40,388) = - reac_rate_local(388) 
  reac_source_local(42,388) = + reac_rate_local(388) 
  reac_source_local(20,389) = - reac_rate_local(389) 
  reac_source_local(40,389) = - reac_rate_local(389) 
  reac_source_local(42,389) = + reac_rate_local(389) 
  reac_source_local(01,390) = + reac_rate_local(390) 
  reac_source_local(19,390) = - reac_rate_local(390) 
  reac_source_local(24,390) = + reac_rate_local(390) 
  reac_source_local(39,390) = - reac_rate_local(390) 
  reac_source_local(01,391) = + reac_rate_local(391) 
  reac_source_local(20,391) = - reac_rate_local(391) 
  reac_source_local(24,391) = + reac_rate_local(391) 
  reac_source_local(39,391) = - reac_rate_local(391) 
  reac_source_local(19,392) = - reac_rate_local(392) 
  reac_source_local(39,392) = - reac_rate_local(392) 
  reac_source_local(41,392) = + reac_rate_local(392) 
  reac_source_local(20,393) = - reac_rate_local(393) 
  reac_source_local(39,393) = - reac_rate_local(393) 
  reac_source_local(41,393) = + reac_rate_local(393) 
  reac_source_local(01,394) = + reac_rate_local(394) 
  reac_source_local(15,394) = + reac_rate_local(394) 
  reac_source_local(40,394) = - reac_rate_local(394) 
  reac_source_local(01,395) = + reac_rate_local(395) 
  reac_source_local(15,395) = + reac_rate_local(395) 
  reac_source_local(40,395) = - reac_rate_local(395) 
  reac_source_local(01,396) = + reac_rate_local(396) 
  reac_source_local(24,396) = + reac_rate_local(396) 
  reac_source_local(41,396) = - reac_rate_local(396) 
  reac_source_local(01,397) = + reac_rate_local(397) 
  reac_source_local(24,397) = + reac_rate_local(397) 
  reac_source_local(41,397) = - reac_rate_local(397) 
  reac_source_local(09,398) = + reac_rate_local(398) 
  reac_source_local(22,398) = - reac_rate_local(398) 
  reac_source_local(28,398) = - reac_rate_local(398) 
  reac_source_local(32,398) = + reac_rate_local(398) 
  reac_source_local(22,399) = - reac_rate_local(399) 
  reac_source_local(25,399) = + reac_rate_local(399) 
  reac_source_local(28,399) = + reac_rate_local(399) 
  reac_source_local(36,399) = - reac_rate_local(399) 
  reac_source_local(15,400) = + reac_rate_local(400) 
  reac_source_local(22,400) = - reac_rate_local(400) 
  reac_source_local(32,400) = + reac_rate_local(400) 
  reac_source_local(36,400) = - reac_rate_local(400) 
  reac_source_local(22,401) = - reac_rate_local(401) 
  reac_source_local(23,401) = + reac_rate_local(401) 
  reac_source_local(32,401) = + reac_rate_local(401) 
  reac_source_local(33,401) = - reac_rate_local(401) 
  reac_source_local(09,402) = + reac_rate_local(402) 
  reac_source_local(25,402) = - reac_rate_local(402) 
  reac_source_local(26,402) = - reac_rate_local(402) 
  reac_source_local(36,402) = + reac_rate_local(402) 
  reac_source_local(02,403) = + reac_rate_local(403) 
  reac_source_local(25,403) = - reac_rate_local(403) 
  reac_source_local(26,403) = - reac_rate_local(403) 
  reac_source_local(37,403) = + reac_rate_local(403) 
  reac_source_local(02,404) = + reac_rate_local(404) * 2.d0
  reac_source_local(25,404) = + reac_rate_local(404) 
  reac_source_local(26,404) = - reac_rate_local(404) 
  reac_source_local(34,404) = - reac_rate_local(404) 
  reac_source_local(02,405) = + reac_rate_local(405) * 2.d0
  reac_source_local(15,405) = + reac_rate_local(405) 
  reac_source_local(26,405) = - reac_rate_local(405) 
  reac_source_local(35,405) = - reac_rate_local(405) 
  reac_source_local(09,406) = + reac_rate_local(406) 
  reac_source_local(15,406) = + reac_rate_local(406) 
  reac_source_local(25,406) = + reac_rate_local(406) 
  reac_source_local(26,406) = - reac_rate_local(406) 
  reac_source_local(40,406) = - reac_rate_local(406) 
  reac_source_local(09,407) = + reac_rate_local(407) 
  reac_source_local(25,407) = - reac_rate_local(407) 
  reac_source_local(28,407) = - reac_rate_local(407) 
  reac_source_local(37,407) = + reac_rate_local(407) 
  reac_source_local(09,408) = + reac_rate_local(408) 
  reac_source_local(28,408) = + reac_rate_local(408) 
  reac_source_local(29,408) = - reac_rate_local(408) 
  reac_source_local(02,409) = + reac_rate_local(409) 
  reac_source_local(09,409) = + reac_rate_local(409) * 2.d0
  reac_source_local(25,409) = + reac_rate_local(409) 
  reac_source_local(29,409) = - reac_rate_local(409) 
  reac_source_local(34,409) = - reac_rate_local(409) 
  reac_source_local(02,410) = + reac_rate_local(410) 
  reac_source_local(09,410) = + reac_rate_local(410) * 2.d0
  reac_source_local(15,410) = + reac_rate_local(410) 
  reac_source_local(29,410) = - reac_rate_local(410) 
  reac_source_local(35,410) = - reac_rate_local(410) 
  reac_source_local(09,411) = + reac_rate_local(411) * 2.d0
  reac_source_local(15,411) = + reac_rate_local(411) 
  reac_source_local(29,411) = - reac_rate_local(411) 
  reac_source_local(40,411) = - reac_rate_local(411) 
  reac_source_local(02,412) = + reac_rate_local(412) * 2.d0
  reac_source_local(09,412) = + reac_rate_local(412) 
  reac_source_local(25,412) = + reac_rate_local(412) 
  reac_source_local(30,412) = - reac_rate_local(412) 
  reac_source_local(34,412) = - reac_rate_local(412) 
  reac_source_local(02,413) = + reac_rate_local(413) * 2.d0
  reac_source_local(09,413) = + reac_rate_local(413) 
  reac_source_local(15,413) = + reac_rate_local(413) 
  reac_source_local(30,413) = - reac_rate_local(413) 
  reac_source_local(35,413) = - reac_rate_local(413) 
  reac_source_local(02,414) = + reac_rate_local(414) 
  reac_source_local(09,414) = + reac_rate_local(414) 
  reac_source_local(15,414) = + reac_rate_local(414) 
  reac_source_local(30,414) = - reac_rate_local(414) 
  reac_source_local(40,414) = - reac_rate_local(414) 
  reac_source_local(02,415) = + reac_rate_local(415) 
  reac_source_local(26,415) = + reac_rate_local(415) 
  reac_source_local(31,415) = - reac_rate_local(415) 
  reac_source_local(02,416) = + reac_rate_local(416) * 3.d0
  reac_source_local(25,416) = + reac_rate_local(416) 
  reac_source_local(31,416) = - reac_rate_local(416) 
  reac_source_local(34,416) = - reac_rate_local(416) 
  reac_source_local(02,417) = + reac_rate_local(417) * 3.d0
  reac_source_local(15,417) = + reac_rate_local(417) 
  reac_source_local(31,417) = - reac_rate_local(417) 
  reac_source_local(35,417) = - reac_rate_local(417) 
  reac_source_local(02,418) = + reac_rate_local(418) * 2.d0
  reac_source_local(15,418) = + reac_rate_local(418) 
  reac_source_local(31,418) = - reac_rate_local(418) 
  reac_source_local(40,418) = - reac_rate_local(418) 
  reac_source_local(02,419) = + reac_rate_local(419) 
  reac_source_local(15,419) = + reac_rate_local(419) 
  reac_source_local(25,419) = + reac_rate_local(419) 
  reac_source_local(34,419) = - reac_rate_local(419) 
  reac_source_local(36,419) = - reac_rate_local(419) 
  reac_source_local(02,420) = + reac_rate_local(420) 
  reac_source_local(25,420) = - reac_rate_local(420) 
  reac_source_local(34,420) = - reac_rate_local(420) 
  reac_source_local(40,420) = + reac_rate_local(420) 
  reac_source_local(02,421) = + reac_rate_local(421) 
  reac_source_local(15,421) = + reac_rate_local(421) * 2.d0
  reac_source_local(35,421) = - reac_rate_local(421) 
  reac_source_local(36,421) = - reac_rate_local(421) 
  reac_source_local(15,422) = + reac_rate_local(422) 
  reac_source_local(25,422) = - reac_rate_local(422) 
  reac_source_local(34,422) = + reac_rate_local(422) 
  reac_source_local(35,422) = - reac_rate_local(422) 
  reac_source_local(02,423) = + reac_rate_local(423) 
  reac_source_local(15,423) = + reac_rate_local(423) 
  reac_source_local(25,423) = - reac_rate_local(423) 
  reac_source_local(35,423) = - reac_rate_local(423) 
  reac_source_local(39,423) = + reac_rate_local(423) 
  reac_source_local(02,424) = + reac_rate_local(424) 
  reac_source_local(25,424) = - reac_rate_local(424) 
  reac_source_local(35,424) = - reac_rate_local(424) 
  reac_source_local(41,424) = + reac_rate_local(424) 
  reac_source_local(02,425) = + reac_rate_local(425) 
  reac_source_local(15,425) = + reac_rate_local(425) 
  reac_source_local(24,425) = - reac_rate_local(425) 
  reac_source_local(35,425) = - reac_rate_local(425) 
  reac_source_local(41,425) = + reac_rate_local(425) 
  reac_source_local(24,426) = - reac_rate_local(426) 
  reac_source_local(25,426) = + reac_rate_local(426) 
  reac_source_local(39,426) = - reac_rate_local(426) 
  reac_source_local(41,426) = + reac_rate_local(426) 
  reac_source_local(01,427) = + reac_rate_local(427) 
  reac_source_local(15,427) = + reac_rate_local(427) * 2.d0
  reac_source_local(24,427) = - reac_rate_local(427) 
  reac_source_local(39,427) = - reac_rate_local(427) 
  reac_source_local(15,428) = + reac_rate_local(428) 
  reac_source_local(24,428) = - reac_rate_local(428) 
  reac_source_local(40,428) = - reac_rate_local(428) 
  reac_source_local(41,428) = + reac_rate_local(428) 
  reac_source_local(01,429) = + reac_rate_local(429) 
  reac_source_local(15,429) = + reac_rate_local(429) * 3.d0
  reac_source_local(24,429) = - reac_rate_local(429) 
  reac_source_local(41,429) = - reac_rate_local(429) 
  reac_source_local(15,430) = + reac_rate_local(430) 
  reac_source_local(24,430) = - reac_rate_local(430) 
  reac_source_local(36,430) = + reac_rate_local(430) 
  reac_source_local(37,430) = - reac_rate_local(430) 
  reac_source_local(25,431) = - reac_rate_local(431) 
  reac_source_local(36,431) = + reac_rate_local(431) 
  reac_source_local(37,431) = - reac_rate_local(431) 
  reac_source_local(01,432) = + reac_rate_local(432) 
  reac_source_local(15,432) = + reac_rate_local(432) 
  reac_source_local(25,432) = - reac_rate_local(432) 
  reac_source_local(39,432) = - reac_rate_local(432) 
  reac_source_local(15,433) = + reac_rate_local(433) 
  reac_source_local(25,433) = - reac_rate_local(433) 
  reac_source_local(39,433) = + reac_rate_local(433) 
  reac_source_local(40,433) = - reac_rate_local(433) 
  reac_source_local(01,434) = + reac_rate_local(434) 
  reac_source_local(24,434) = + reac_rate_local(434) 
  reac_source_local(25,434) = - reac_rate_local(434) 
  reac_source_local(40,434) = - reac_rate_local(434) 
  reac_source_local(24,435) = + reac_rate_local(435) 
  reac_source_local(25,435) = - reac_rate_local(435) 
  reac_source_local(39,435) = + reac_rate_local(435) 
  reac_source_local(41,435) = - reac_rate_local(435) 
  reac_source_local(01,436) = + reac_rate_local(436) 
  reac_source_local(15,436) = + reac_rate_local(436) * 2.d0
  reac_source_local(25,436) = - reac_rate_local(436) 
  reac_source_local(41,436) = - reac_rate_local(436) 
  reac_source_local(15,437) = + reac_rate_local(437) 
  reac_source_local(25,437) = - reac_rate_local(437) 
  reac_source_local(40,437) = + reac_rate_local(437) 
  reac_source_local(41,437) = - reac_rate_local(437) 
  reac_source_local(15,438) = + reac_rate_local(438) 
  reac_source_local(25,438) = - reac_rate_local(438) 
  reac_source_local(41,438) = + reac_rate_local(438) 
  reac_source_local(42,438) = - reac_rate_local(438) 
  reac_source_local(15,439) = + reac_rate_local(439) * 2.d0
  reac_source_local(25,439) = - reac_rate_local(439) 
  reac_source_local(39,439) = + reac_rate_local(439) 
  reac_source_local(42,439) = - reac_rate_local(439) 
  reac_source_local(24,440) = + reac_rate_local(440) 
  reac_source_local(25,440) = - reac_rate_local(440) 
  reac_source_local(40,440) = + reac_rate_local(440) 
  reac_source_local(42,440) = - reac_rate_local(440) 
  reac_source_local(24,441) = + reac_rate_local(441) 
  reac_source_local(37,441) = - reac_rate_local(441) 
  reac_source_local(40,441) = - reac_rate_local(441) 
  reac_source_local(15,442) = + reac_rate_local(442) 
  reac_source_local(25,442) = + reac_rate_local(442) 
  reac_source_local(37,442) = - reac_rate_local(442) 
  reac_source_local(40,442) = - reac_rate_local(442) 
  reac_source_local(15,443) = + reac_rate_local(443) * 2.d0
  reac_source_local(36,443) = - reac_rate_local(443) 
  reac_source_local(40,443) = - reac_rate_local(443) 
  reac_source_local(15,444) = + reac_rate_local(444) 
  reac_source_local(25,444) = + reac_rate_local(444) * 2.d0
  reac_source_local(36,444) = - reac_rate_local(444) 
  reac_source_local(40,444) = - reac_rate_local(444) 
  reac_source_local(15,445) = + reac_rate_local(445) * 2.d0
  reac_source_local(36,445) = - reac_rate_local(445) 
  reac_source_local(40,445) = - reac_rate_local(445) 
  reac_source_local(01,446) = + reac_rate_local(446) 
  reac_source_local(15,446) = + reac_rate_local(446) 
  reac_source_local(40,446) = - reac_rate_local(446) 
  reac_source_local(15,447) = + reac_rate_local(447) 
  reac_source_local(24,447) = + reac_rate_local(447) 
  reac_source_local(36,447) = - reac_rate_local(447) 
  reac_source_local(41,447) = - reac_rate_local(447) 
  reac_source_local(24,448) = + reac_rate_local(448) 
  reac_source_local(25,448) = + reac_rate_local(448) * 2.d0
  reac_source_local(36,448) = - reac_rate_local(448) 
  reac_source_local(41,448) = - reac_rate_local(448) 
  reac_source_local(24,449) = + reac_rate_local(449) 
  reac_source_local(25,449) = + reac_rate_local(449) 
  reac_source_local(37,449) = - reac_rate_local(449) 
  reac_source_local(41,449) = - reac_rate_local(449) 
  reac_source_local(01,450) = + reac_rate_local(450) 
  reac_source_local(24,450) = + reac_rate_local(450) 
  reac_source_local(41,450) = - reac_rate_local(450) 
  reac_source_local(25,451) = + reac_rate_local(451) * 2.d0
  reac_source_local(37,451) = - reac_rate_local(451) 
  reac_source_local(39,451) = - reac_rate_local(451) 
  reac_source_local(15,452) = + reac_rate_local(452) 
  reac_source_local(37,452) = - reac_rate_local(452) 
  reac_source_local(39,452) = - reac_rate_local(452) 
  reac_source_local(15,453) = + reac_rate_local(453) 
  reac_source_local(25,453) = + reac_rate_local(453) 
  reac_source_local(36,453) = - reac_rate_local(453) 
  reac_source_local(39,453) = - reac_rate_local(453) 
  reac_source_local(25,454) = + reac_rate_local(454) * 3.d0
  reac_source_local(36,454) = - reac_rate_local(454) 
  reac_source_local(39,454) = - reac_rate_local(454) 
  reac_source_local(24,455) = + reac_rate_local(455) 
  reac_source_local(36,455) = - reac_rate_local(455) 
  reac_source_local(39,455) = - reac_rate_local(455) 
  reac_source_local(01,456) = + reac_rate_local(456) 
  reac_source_local(25,456) = + reac_rate_local(456) 
  reac_source_local(39,456) = - reac_rate_local(456) 
  reac_source_local(15,457) = + reac_rate_local(457) 
  reac_source_local(40,457) = + reac_rate_local(457) 
  reac_source_local(42,457) = - reac_rate_local(457) 
  reac_source_local(15,458) = + reac_rate_local(458) 
  reac_source_local(36,458) = + reac_rate_local(458) 
  reac_source_local(38,458) = - reac_rate_local(458) 
  reac_source_local(02,459) = + reac_rate_local(459) 
  reac_source_local(03,459) = - reac_rate_local(459) 
  reac_source_local(02,460) = + reac_rate_local(460) 
  reac_source_local(03,460) = - reac_rate_local(460) 
  reac_source_local(02,461) = + reac_rate_local(461) 
  reac_source_local(03,461) = - reac_rate_local(461) 
  reac_source_local(02,462) = + reac_rate_local(462) 
  reac_source_local(04,462) = - reac_rate_local(462) 
  reac_source_local(02,463) = + reac_rate_local(463) 
  reac_source_local(04,463) = - reac_rate_local(463) 
  reac_source_local(02,464) = + reac_rate_local(464) 
  reac_source_local(04,464) = - reac_rate_local(464) 
  reac_source_local(03,465) = + reac_rate_local(465) 
  reac_source_local(04,465) = - reac_rate_local(465) 
  reac_source_local(03,466) = + reac_rate_local(466) 
  reac_source_local(04,466) = - reac_rate_local(466) 
  reac_source_local(03,467) = + reac_rate_local(467) 
  reac_source_local(04,467) = - reac_rate_local(467) 
  reac_source_local(04,468) = + reac_rate_local(468) 
  reac_source_local(05,468) = - reac_rate_local(468) 
  reac_source_local(04,469) = + reac_rate_local(469) 
  reac_source_local(05,469) = - reac_rate_local(469) 
  reac_source_local(04,470) = + reac_rate_local(470) 
  reac_source_local(05,470) = - reac_rate_local(470) 
  reac_source_local(05,471) = - reac_rate_local(471) 
  reac_source_local(06,471) = + reac_rate_local(471) 
  reac_source_local(05,472) = - reac_rate_local(472) 
  reac_source_local(06,472) = + reac_rate_local(472) 
  reac_source_local(05,473) = - reac_rate_local(473) 
  reac_source_local(06,473) = + reac_rate_local(473) 
  reac_source_local(02,474) = - reac_rate_local(474) 
  reac_source_local(03,474) = + reac_rate_local(474) 
  reac_source_local(04,474) = + reac_rate_local(474) 
  reac_source_local(05,474) = - reac_rate_local(474) 
  reac_source_local(03,475) = + reac_rate_local(475) 
  reac_source_local(05,475) = - reac_rate_local(475) 
  reac_source_local(03,476) = + reac_rate_local(476) 
  reac_source_local(05,476) = - reac_rate_local(476) 
  reac_source_local(03,477) = + reac_rate_local(477) 
  reac_source_local(05,477) = - reac_rate_local(477) 
  reac_source_local(04,478) = + reac_rate_local(478) 
  reac_source_local(06,478) = - reac_rate_local(478) 
  reac_source_local(04,479) = + reac_rate_local(479) 
  reac_source_local(06,479) = - reac_rate_local(479) 
  reac_source_local(04,480) = + reac_rate_local(480) 
  reac_source_local(06,480) = - reac_rate_local(480) 
  reac_source_local(03,481) = + reac_rate_local(481) 
  reac_source_local(06,481) = - reac_rate_local(481) 
  reac_source_local(03,482) = + reac_rate_local(482) 
  reac_source_local(06,482) = - reac_rate_local(482) 
  reac_source_local(03,483) = + reac_rate_local(483) 
  reac_source_local(06,483) = - reac_rate_local(483) 
  reac_source_local(09,484) = + reac_rate_local(484) 
  reac_source_local(10,484) = - reac_rate_local(484) 
  reac_source_local(09,485) = + reac_rate_local(485) 
  reac_source_local(10,485) = - reac_rate_local(485) 
  reac_source_local(09,486) = + reac_rate_local(486) 
  reac_source_local(10,486) = - reac_rate_local(486) 
  reac_source_local(15,487) = + reac_rate_local(487) 
  reac_source_local(16,487) = - reac_rate_local(487) 
  reac_source_local(15,488) = + reac_rate_local(488) 
  reac_source_local(16,488) = - reac_rate_local(488) 
  reac_source_local(15,489) = + reac_rate_local(489) 
  reac_source_local(16,489) = - reac_rate_local(489) 
  reac_source_local(15,490) = + reac_rate_local(490) 
  reac_source_local(17,490) = - reac_rate_local(490) 
  reac_source_local(15,491) = + reac_rate_local(491) 
  reac_source_local(17,491) = - reac_rate_local(491) 
  reac_source_local(15,492) = + reac_rate_local(492) 
  reac_source_local(17,492) = - reac_rate_local(492) 
  reac_source_local(15,493) = + reac_rate_local(493) 
  reac_source_local(18,493) = - reac_rate_local(493) 
  reac_source_local(15,494) = + reac_rate_local(494) 
  reac_source_local(18,494) = - reac_rate_local(494) 
  reac_source_local(15,495) = + reac_rate_local(495) 
  reac_source_local(18,495) = - reac_rate_local(495) 
  reac_source_local(02,496) = + reac_rate_local(496) 
  reac_source_local(07,496) = - reac_rate_local(496) 
  reac_source_local(02,497) = + reac_rate_local(497) 
  reac_source_local(08,497) = - reac_rate_local(497) 
  reac_source_local(09,498) = + reac_rate_local(498) 
  reac_source_local(11,498) = - reac_rate_local(498) 
  reac_source_local(09,499) = + reac_rate_local(499) 
  reac_source_local(12,499) = - reac_rate_local(499) 
  reac_source_local(09,500) = + reac_rate_local(500) 
  reac_source_local(13,500) = - reac_rate_local(500) 
  reac_source_local(09,501) = + reac_rate_local(501) 
  reac_source_local(14,501) = - reac_rate_local(501) 
  reac_source_local(15,502) = + reac_rate_local(502) 
  reac_source_local(19,502) = - reac_rate_local(502) 
  reac_source_local(15,503) = + reac_rate_local(503) 
  reac_source_local(20,503) = - reac_rate_local(503) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(43)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(02) 
  rrt(002) = rrt(002) * density(01) * density(02) 
  rrt(003) = rrt(003) * density(01) * density(02) 
  rrt(004) = rrt(004) * density(01) * density(02) 
  rrt(005) = rrt(005) * density(01) * density(02) 
  rrt(006) = rrt(006) * density(01) * density(02) 
  rrt(007) = rrt(007) * density(01) * density(02) 
  rrt(008) = rrt(008) * density(01) * density(02) 
  rrt(009) = rrt(009) * density(01) * density(02) 
  rrt(010) = rrt(010) * density(01) * density(02) 
  rrt(011) = rrt(011) * density(01) * density(02) 
  rrt(012) = rrt(012) * density(01) * density(02) 
  rrt(013) = rrt(013) * density(01) * density(02) 
  rrt(014) = rrt(014) * density(01) * density(07) 
  rrt(015) = rrt(015) * density(01) * density(08) 
  rrt(016) = rrt(016) * density(01) * density(03) 
  rrt(017) = rrt(017) * density(01) * density(04) 
  rrt(018) = rrt(018) * density(01) * density(06) 
  rrt(019) = rrt(019) * density(01) * density(05) 
  rrt(020) = rrt(020) * density(01) * density(03) 
  rrt(021) = rrt(021) * density(01) * density(04) 
  rrt(022) = rrt(022) * density(01) * density(06) 
  rrt(023) = rrt(023) * density(01) * density(05) 
  rrt(024) = rrt(024) * density(01) * density(03) 
  rrt(025) = rrt(025) * density(01) * density(04) 
  rrt(026) = rrt(026) * density(01) * density(06) 
  rrt(027) = rrt(027) * density(01) * density(05) 
  rrt(028) = rrt(028) * density(01) * density(03) 
  rrt(029) = rrt(029) * density(01) * density(04) 
  rrt(030) = rrt(030) * density(01) * density(06) 
  rrt(031) = rrt(031) * density(01) * density(05) 
  rrt(032) = rrt(032) * density(01) * density(03) 
  rrt(033) = rrt(033) * density(01) * density(04) 
  rrt(034) = rrt(034) * density(01) * density(06) 
  rrt(035) = rrt(035) * density(01) * density(05) 
  rrt(036) = rrt(036) * density(01) * density(03) 
  rrt(037) = rrt(037) * density(01) * density(04) 
  rrt(038) = rrt(038) * density(01) * density(06) 
  rrt(039) = rrt(039) * density(01) * density(05) 
  rrt(040) = rrt(040) * density(01) * density(03) 
  rrt(041) = rrt(041) * density(01) * density(04) 
  rrt(042) = rrt(042) * density(01) * density(06) 
  rrt(043) = rrt(043) * density(01) * density(05) 
  rrt(044) = rrt(044) * density(01) * density(03) 
  rrt(045) = rrt(045) * density(01) * density(04) 
  rrt(046) = rrt(046) * density(01) * density(06) 
  rrt(047) = rrt(047) * density(01) * density(05) 
  rrt(048) = rrt(048) * density(01) * density(03) 
  rrt(049) = rrt(049) * density(01) * density(04) 
  rrt(050) = rrt(050) * density(01) * density(06) 
  rrt(051) = rrt(051) * density(01) * density(05) 
  rrt(052) = rrt(052) * density(01) * density(03) 
  rrt(053) = rrt(053) * density(01) * density(04) 
  rrt(054) = rrt(054) * density(01) * density(06) 
  rrt(055) = rrt(055) * density(01) * density(05) 
  rrt(056) = rrt(056) * density(01) * density(07) 
  rrt(057) = rrt(057) * density(01) * density(08) 
  rrt(058) = rrt(058) * density(01) * density(07) 
  rrt(059) = rrt(059) * density(01) * density(08) 
  rrt(060) = rrt(060) * density(01) * density(07) 
  rrt(061) = rrt(061) * density(01) * density(08) 
  rrt(062) = rrt(062) * density(01) * density(07) 
  rrt(063) = rrt(063) * density(01) * density(08) 
  rrt(064) = rrt(064) * density(01) * density(07) 
  rrt(065) = rrt(065) * density(01) * density(08) 
  rrt(066) = rrt(066) * density(01) * density(07) 
  rrt(067) = rrt(067) * density(01) * density(08) 
  rrt(068) = rrt(068) * density(01) * density(09) 
  rrt(069) = rrt(069) * density(01) * density(09) 
  rrt(070) = rrt(070) * density(01) * density(09) 
  rrt(071) = rrt(071) * density(01) * density(09) 
  rrt(072) = rrt(072) * density(01) * density(09) 
  rrt(073) = rrt(073) * density(01) * density(09) 
  rrt(074) = rrt(074) * density(01) * density(09) 
  rrt(075) = rrt(075) * density(01) * density(09) 
  rrt(076) = rrt(076) * density(01) * density(09) 
  rrt(077) = rrt(077) * density(01) * density(09) 
  rrt(078) = rrt(078) * density(01) * density(11) 
  rrt(079) = rrt(079) * density(01) * density(12) 
  rrt(080) = rrt(080) * density(01) * density(13) 
  rrt(081) = rrt(081) * density(01) * density(14) 
  rrt(082) = rrt(082) * density(01) * density(10) 
  rrt(083) = rrt(083) * density(01) * density(10) 
  rrt(084) = rrt(084) * density(01) * density(10) 
  rrt(085) = rrt(085) * density(01) * density(10) 
  rrt(086) = rrt(086) * density(01) * density(10) 
  rrt(087) = rrt(087) * density(01) * density(10) 
  rrt(088) = rrt(088) * density(01) * density(10) 
  rrt(089) = rrt(089) * density(01) * density(10) 
  rrt(090) = rrt(090) * density(01) * density(10) 
  rrt(091) = rrt(091) * density(01) * density(10) 
  rrt(092) = rrt(092) * density(01) * density(11) 
  rrt(093) = rrt(093) * density(01) * density(12) 
  rrt(094) = rrt(094) * density(01) * density(13) 
  rrt(095) = rrt(095) * density(01) * density(14) 
  rrt(096) = rrt(096) * density(01) * density(11) 
  rrt(097) = rrt(097) * density(01) * density(12) 
  rrt(098) = rrt(098) * density(01) * density(13) 
  rrt(099) = rrt(099) * density(01) * density(14) 
  rrt(100) = rrt(100) * density(01) * density(11) 
  rrt(101) = rrt(101) * density(01) * density(12) 
  rrt(102) = rrt(102) * density(01) * density(13) 
  rrt(103) = rrt(103) * density(01) * density(14) 
  rrt(104) = rrt(104) * density(01) * density(11) 
  rrt(105) = rrt(105) * density(01) * density(12) 
  rrt(106) = rrt(106) * density(01) * density(13) 
  rrt(107) = rrt(107) * density(01) * density(14) 
  rrt(108) = rrt(108) * density(01) * density(22) 
  rrt(109) = rrt(109) * density(01) * density(23) 
  rrt(110) = rrt(110) * density(01) * density(23) 
  rrt(111) = rrt(111) * density(01) * density(15) 
  rrt(112) = rrt(112) * density(01) * density(15) 
  rrt(113) = rrt(113) * density(01) * density(15) 
  rrt(114) = rrt(114) * density(01) * density(15) 
  rrt(115) = rrt(115) * density(01) * density(15) 
  rrt(116) = rrt(116) * density(01) * density(15) 
  rrt(117) = rrt(117) * density(01) * density(15) 
  rrt(118) = rrt(118) * density(01) * density(15) 
  rrt(119) = rrt(119) * density(01) * density(15) 
  rrt(120) = rrt(120) * density(01) * density(19) 
  rrt(121) = rrt(121) * density(01) * density(20) 
  rrt(122) = rrt(122) * density(19) 
  rrt(123) = rrt(123) * density(20) 
  rrt(124) = rrt(124) * density(01) * density(16) 
  rrt(125) = rrt(125) * density(01) * density(17) 
  rrt(126) = rrt(126) * density(01) * density(18) 
  rrt(127) = rrt(127) * density(01) * density(16) 
  rrt(128) = rrt(128) * density(01) * density(17) 
  rrt(129) = rrt(129) * density(01) * density(18) 
  rrt(130) = rrt(130) * density(01) * density(16) 
  rrt(131) = rrt(131) * density(01) * density(17) 
  rrt(132) = rrt(132) * density(01) * density(18) 
  rrt(133) = rrt(133) * density(01) * density(16) 
  rrt(134) = rrt(134) * density(01) * density(17) 
  rrt(135) = rrt(135) * density(01) * density(18) 
  rrt(136) = rrt(136) * density(01) * density(16) 
  rrt(137) = rrt(137) * density(01) * density(17) 
  rrt(138) = rrt(138) * density(01) * density(18) 
  rrt(139) = rrt(139) * density(01) * density(16) 
  rrt(140) = rrt(140) * density(01) * density(17) 
  rrt(141) = rrt(141) * density(01) * density(18) 
  rrt(142) = rrt(142) * density(01) * density(16) 
  rrt(143) = rrt(143) * density(01) * density(17) 
  rrt(144) = rrt(144) * density(01) * density(18) 
  rrt(145) = rrt(145) * density(01) * density(19) 
  rrt(146) = rrt(146) * density(01) * density(20) 
  rrt(147) = rrt(147) * density(01) * density(19) 
  rrt(148) = rrt(148) * density(01) * density(20) 
  rrt(149) = rrt(149) * density(01) * density(19) 
  rrt(150) = rrt(150) * density(01) * density(20) 
  rrt(151) = rrt(151) * density(01) * density(24) 
  rrt(152) = rrt(152) * density(01) * density(24) 
  rrt(153) = rrt(153) * density(01) * density(24) 
  rrt(154) = rrt(154) * density(01) * density(24) 
  rrt(155) = rrt(155) * density(01) * density(24) 
  rrt(156) = rrt(156) * density(01) * density(25) 
  rrt(157) = rrt(157) * density(01) * density(26) 
  rrt(158) = rrt(158) * density(01) * density(26) 
  rrt(159) = rrt(159) * density(01) * density(27) 
  rrt(160) = rrt(160) * density(01) * density(28) 
  rrt(161) = rrt(161) * density(01) * density(29) 
  rrt(162) = rrt(162) * density(01) * density(30) 
  rrt(163) = rrt(163) * density(01) * density(31) 
  rrt(164) = rrt(164) * density(01) * density(33) 
  rrt(165) = rrt(165) * density(01) * density(15) 
  rrt(166) = rrt(166) * density(01) * density(24) 
  rrt(167) = rrt(167) * density(01) * density(25) 
  rrt(168) = rrt(168) * density(01) * density(36) 
  rrt(169) = rrt(169) * density(01) * density(36) 
  rrt(170) = rrt(170) * density(01) * density(37) 
  rrt(171) = rrt(171) * density(01) * density(38) 
  rrt(172) = rrt(172) * density(02) * density(25) 
  rrt(173) = rrt(173) * density(02) * density(22) 
  rrt(174) = rrt(174) * density(03) * density(25) 
  rrt(175) = rrt(175) * density(04) * density(25) 
  rrt(176) = rrt(176) * density(05) * density(25) 
  rrt(177) = rrt(177) * density(06) * density(25) 
  rrt(178) = rrt(178) * density(03) * density(22) 
  rrt(179) = rrt(179) * density(04) * density(22) 
  rrt(180) = rrt(180) * density(05) * density(22) 
  rrt(181) = rrt(181) * density(06) * density(22) 
  rrt(182) = rrt(182) * density(09) * density(25) 
  rrt(183) = rrt(183) * density(09) * density(24) 
  rrt(184) = rrt(184) * density(09) * density(22) 
  rrt(185) = rrt(185) * density(09) * density(15) 
  rrt(186) = rrt(186) * density(10) * density(25) 
  rrt(187) = rrt(187) * density(09) * density(16) 
  rrt(188) = rrt(188) * density(09) * density(17) 
  rrt(189) = rrt(189) * density(09) * density(18) 
  rrt(190) = rrt(190) * density(10) * density(15) 
  rrt(191) = rrt(191) * density(10) * density(16) 
  rrt(192) = rrt(192) * density(10) * density(17) 
  rrt(193) = rrt(193) * density(10) * density(18) 
  rrt(194) = rrt(194) * density(10) * density(24) 
  rrt(195) = rrt(195) * density(10) * density(22) 
  rrt(196) = rrt(196) * density(15) * density(22) 
  rrt(197) = rrt(197) * density(15) * density(21) 
  rrt(198) = rrt(198) * density(15) * density(25) 
  rrt(199) = rrt(199) * density(16) * density(22) 
  rrt(200) = rrt(200) * density(17) * density(22) 
  rrt(201) = rrt(201) * density(18) * density(22) 
  rrt(202) = rrt(202) * density(16) * density(25) 
  rrt(203) = rrt(203) * density(17) * density(25) 
  rrt(204) = rrt(204) * density(18) * density(25) 
  rrt(205) = rrt(205) * density(22) * density(25) 
  rrt(206) = rrt(206) * density(21) * density(25) 
  rrt(207) = rrt(207) * density(24) * density(25) 
  rrt(208) = rrt(208) * density(24) 
  rrt(209) = rrt(209) * density(25)**2 
  rrt(210) = rrt(210) * density(02) * density(36) 
  rrt(211) = rrt(211) * density(02) * density(37) 
  rrt(212) = rrt(212) * density(02) * density(37) 
  rrt(213) = rrt(213) * density(02) * density(32) 
  rrt(214) = rrt(214) * density(02) * density(28) 
  rrt(215) = rrt(215) * density(02) * density(39) 
  rrt(216) = rrt(216) * density(02) * density(40) 
  rrt(217) = rrt(217) * density(02) * density(41) 
  rrt(218) = rrt(218) * density(02) * density(42) 
  rrt(219) = rrt(219) * density(02) * density(26) 
  rrt(220) = rrt(220) * density(03) * density(37) 
  rrt(221) = rrt(221) * density(04) * density(37) 
  rrt(222) = rrt(222) * density(05) * density(37) 
  rrt(223) = rrt(223) * density(06) * density(37) 
  rrt(224) = rrt(224) * density(03) * density(37) 
  rrt(225) = rrt(225) * density(04) * density(37) 
  rrt(226) = rrt(226) * density(05) * density(37) 
  rrt(227) = rrt(227) * density(06) * density(37) 
  rrt(228) = rrt(228) * density(03) * density(32) 
  rrt(229) = rrt(229) * density(04) * density(32) 
  rrt(230) = rrt(230) * density(05) * density(32) 
  rrt(231) = rrt(231) * density(06) * density(32) 
  rrt(232) = rrt(232) * density(03) * density(28) 
  rrt(233) = rrt(233) * density(04) * density(28) 
  rrt(234) = rrt(234) * density(05) * density(28) 
  rrt(235) = rrt(235) * density(06) * density(28) 
  rrt(236) = rrt(236) * density(03) * density(39) 
  rrt(237) = rrt(237) * density(04) * density(39) 
  rrt(238) = rrt(238) * density(05) * density(39) 
  rrt(239) = rrt(239) * density(06) * density(39) 
  rrt(240) = rrt(240) * density(03) * density(40) 
  rrt(241) = rrt(241) * density(04) * density(40) 
  rrt(242) = rrt(242) * density(05) * density(40) 
  rrt(243) = rrt(243) * density(06) * density(40) 
  rrt(244) = rrt(244) * density(03) * density(41) 
  rrt(245) = rrt(245) * density(04) * density(41) 
  rrt(246) = rrt(246) * density(05) * density(41) 
  rrt(247) = rrt(247) * density(06) * density(41) 
  rrt(248) = rrt(248) * density(03) * density(42) 
  rrt(249) = rrt(249) * density(04) * density(42) 
  rrt(250) = rrt(250) * density(05) * density(42) 
  rrt(251) = rrt(251) * density(06) * density(42) 
  rrt(252) = rrt(252) * density(03) * density(26) 
  rrt(253) = rrt(253) * density(04) * density(26) 
  rrt(254) = rrt(254) * density(05) * density(26) 
  rrt(255) = rrt(255) * density(06) * density(26) 
  rrt(256) = rrt(256) * density(03) * density(36) 
  rrt(257) = rrt(257) * density(04) * density(36) 
  rrt(258) = rrt(258) * density(05) * density(36) 
  rrt(259) = rrt(259) * density(06) * density(36) 
  rrt(260) = rrt(260) * density(07) * density(37) 
  rrt(261) = rrt(261) * density(08) * density(37) 
  rrt(262) = rrt(262) * density(07) * density(37) 
  rrt(263) = rrt(263) * density(08) * density(37) 
  rrt(264) = rrt(264) * density(07) * density(32) 
  rrt(265) = rrt(265) * density(08) * density(32) 
  rrt(266) = rrt(266) * density(07) * density(28) 
  rrt(267) = rrt(267) * density(08) * density(28) 
  rrt(268) = rrt(268) * density(07) * density(39) 
  rrt(269) = rrt(269) * density(08) * density(39) 
  rrt(270) = rrt(270) * density(07) * density(40) 
  rrt(271) = rrt(271) * density(08) * density(40) 
  rrt(272) = rrt(272) * density(07) * density(41) 
  rrt(273) = rrt(273) * density(08) * density(41) 
  rrt(274) = rrt(274) * density(07) * density(42) 
  rrt(275) = rrt(275) * density(08) * density(42) 
  rrt(276) = rrt(276) * density(07) * density(26) 
  rrt(277) = rrt(277) * density(08) * density(26) 
  rrt(278) = rrt(278) * density(07) * density(36) 
  rrt(279) = rrt(279) * density(08) * density(36) 
  rrt(280) = rrt(280) * density(07) * density(36) 
  rrt(281) = rrt(281) * density(08) * density(36) 
  rrt(282) = rrt(282) * density(09) * density(37) 
  rrt(283) = rrt(283) * density(09) * density(39) 
  rrt(284) = rrt(284) * density(09) * density(34) 
  rrt(285) = rrt(285) * density(09) * density(30) 
  rrt(286) = rrt(286) * density(09) * density(31) 
  rrt(287) = rrt(287) * density(09) * density(30) 
  rrt(288) = rrt(288) * density(09) * density(31) 
  rrt(289) = rrt(289) * density(09) * density(32) 
  rrt(290) = rrt(290) * density(10) * density(37) 
  rrt(291) = rrt(291) * density(10) * density(39) 
  rrt(292) = rrt(292) * density(10) * density(34) 
  rrt(293) = rrt(293) * density(10) * density(30) 
  rrt(294) = rrt(294) * density(10) * density(31) 
  rrt(295) = rrt(295) * density(10) * density(30) 
  rrt(296) = rrt(296) * density(10) * density(31) 
  rrt(297) = rrt(297) * density(10) * density(32) 
  rrt(298) = rrt(298) * density(11) * density(39) 
  rrt(299) = rrt(299) * density(12) * density(39) 
  rrt(300) = rrt(300) * density(13) * density(39) 
  rrt(301) = rrt(301) * density(14) * density(39) 
  rrt(302) = rrt(302) * density(11) * density(34) 
  rrt(303) = rrt(303) * density(12) * density(34) 
  rrt(304) = rrt(304) * density(13) * density(34) 
  rrt(305) = rrt(305) * density(14) * density(34) 
  rrt(306) = rrt(306) * density(11) * density(30) 
  rrt(307) = rrt(307) * density(12) * density(30) 
  rrt(308) = rrt(308) * density(13) * density(30) 
  rrt(309) = rrt(309) * density(14) * density(30) 
  rrt(310) = rrt(310) * density(11) * density(31) 
  rrt(311) = rrt(311) * density(12) * density(31) 
  rrt(312) = rrt(312) * density(13) * density(31) 
  rrt(313) = rrt(313) * density(14) * density(31) 
  rrt(314) = rrt(314) * density(11) * density(30) 
  rrt(315) = rrt(315) * density(12) * density(30) 
  rrt(316) = rrt(316) * density(13) * density(30) 
  rrt(317) = rrt(317) * density(14) * density(30) 
  rrt(318) = rrt(318) * density(11) * density(31) 
  rrt(319) = rrt(319) * density(12) * density(31) 
  rrt(320) = rrt(320) * density(13) * density(31) 
  rrt(321) = rrt(321) * density(14) * density(31) 
  rrt(322) = rrt(322) * density(11) * density(32) 
  rrt(323) = rrt(323) * density(12) * density(32) 
  rrt(324) = rrt(324) * density(13) * density(32) 
  rrt(325) = rrt(325) * density(14) * density(32) 
  rrt(326) = rrt(326) * density(15) * density(26) 
  rrt(327) = rrt(327) * density(15) * density(28) 
  rrt(328) = rrt(328) * density(15) * density(29) 
  rrt(329) = rrt(329) * density(15) * density(32) 
  rrt(330) = rrt(330) * density(15) * density(32) 
  rrt(331) = rrt(331) * density(15) * density(37) 
  rrt(332) = rrt(332) * density(15) * density(36) 
  rrt(333) = rrt(333) * density(15) * density(40) 
  rrt(334) = rrt(334) * density(15) * density(39) 
  rrt(335) = rrt(335) * density(15) * density(39) 
  rrt(336) = rrt(336) * density(15) * density(40) 
  rrt(337) = rrt(337) * density(15) * density(41) 
  rrt(338) = rrt(338) * density(16) * density(26) 
  rrt(339) = rrt(339) * density(17) * density(26) 
  rrt(340) = rrt(340) * density(18) * density(26) 
  rrt(341) = rrt(341) * density(16) * density(28) 
  rrt(342) = rrt(342) * density(17) * density(28) 
  rrt(343) = rrt(343) * density(18) * density(28) 
  rrt(344) = rrt(344) * density(16) * density(29) 
  rrt(345) = rrt(345) * density(17) * density(29) 
  rrt(346) = rrt(346) * density(18) * density(29) 
  rrt(347) = rrt(347) * density(16) * density(32) 
  rrt(348) = rrt(348) * density(17) * density(32) 
  rrt(349) = rrt(349) * density(18) * density(32) 
  rrt(350) = rrt(350) * density(16) * density(32) 
  rrt(351) = rrt(351) * density(17) * density(32) 
  rrt(352) = rrt(352) * density(18) * density(32) 
  rrt(353) = rrt(353) * density(16) * density(37) 
  rrt(354) = rrt(354) * density(17) * density(37) 
  rrt(355) = rrt(355) * density(18) * density(37) 
  rrt(356) = rrt(356) * density(16) * density(36) 
  rrt(357) = rrt(357) * density(17) * density(36) 
  rrt(358) = rrt(358) * density(18) * density(36) 
  rrt(359) = rrt(359) * density(16) * density(40) 
  rrt(360) = rrt(360) * density(17) * density(40) 
  rrt(361) = rrt(361) * density(18) * density(40) 
  rrt(362) = rrt(362) * density(16) * density(39) 
  rrt(363) = rrt(363) * density(17) * density(39) 
  rrt(364) = rrt(364) * density(18) * density(39) 
  rrt(365) = rrt(365) * density(16) * density(39) 
  rrt(366) = rrt(366) * density(17) * density(39) 
  rrt(367) = rrt(367) * density(18) * density(39) 
  rrt(368) = rrt(368) * density(16) * density(40) 
  rrt(369) = rrt(369) * density(17) * density(40) 
  rrt(370) = rrt(370) * density(18) * density(40) 
  rrt(371) = rrt(371) * density(16) * density(41) 
  rrt(372) = rrt(372) * density(17) * density(41) 
  rrt(373) = rrt(373) * density(18) * density(41) 
  rrt(374) = rrt(374) * density(19) * density(26) 
  rrt(375) = rrt(375) * density(20) * density(26) 
  rrt(376) = rrt(376) * density(19) * density(28) 
  rrt(377) = rrt(377) * density(20) * density(28) 
  rrt(378) = rrt(378) * density(19) * density(29) 
  rrt(379) = rrt(379) * density(20) * density(29) 
  rrt(380) = rrt(380) * density(19) * density(32) 
  rrt(381) = rrt(381) * density(20) * density(32) 
  rrt(382) = rrt(382) * density(19) * density(32) 
  rrt(383) = rrt(383) * density(20) * density(32) 
  rrt(384) = rrt(384) * density(19) * density(37) 
  rrt(385) = rrt(385) * density(20) * density(37) 
  rrt(386) = rrt(386) * density(19) * density(36) 
  rrt(387) = rrt(387) * density(20) * density(36) 
  rrt(388) = rrt(388) * density(19) * density(40) 
  rrt(389) = rrt(389) * density(20) * density(40) 
  rrt(390) = rrt(390) * density(19) * density(39) 
  rrt(391) = rrt(391) * density(20) * density(39) 
  rrt(392) = rrt(392) * density(19) * density(39) 
  rrt(393) = rrt(393) * density(20) * density(39) 
  rrt(394) = rrt(394) * density(19) * density(40) 
  rrt(395) = rrt(395) * density(20) * density(40) 
  rrt(396) = rrt(396) * density(19) * density(41) 
  rrt(397) = rrt(397) * density(20) * density(41) 
  rrt(398) = rrt(398) * density(22) * density(28) 
  rrt(399) = rrt(399) * density(22) * density(36) 
  rrt(400) = rrt(400) * density(22) * density(36) 
  rrt(401) = rrt(401) * density(22) * density(33) 
  rrt(402) = rrt(402) * density(25) * density(26) 
  rrt(403) = rrt(403) * density(25) * density(26) 
  rrt(404) = rrt(404) * density(26) * density(34) 
  rrt(405) = rrt(405) * density(26) * density(35) 
  rrt(406) = rrt(406) * density(26) * density(40) 
  rrt(407) = rrt(407) * density(25) * density(28) 
  rrt(408) = rrt(408) * density(29) 
  rrt(409) = rrt(409) * density(29) * density(34) 
  rrt(410) = rrt(410) * density(29) * density(35) 
  rrt(411) = rrt(411) * density(29) * density(40) 
  rrt(412) = rrt(412) * density(30) * density(34) 
  rrt(413) = rrt(413) * density(30) * density(35) 
  rrt(414) = rrt(414) * density(30) * density(40) 
  rrt(415) = rrt(415) * density(31) 
  rrt(416) = rrt(416) * density(31) * density(34) 
  rrt(417) = rrt(417) * density(31) * density(35) 
  rrt(418) = rrt(418) * density(31) * density(40) 
  rrt(419) = rrt(419) * density(34) * density(36) 
  rrt(420) = rrt(420) * density(25) * density(34) 
  rrt(421) = rrt(421) * density(35) * density(36) 
  rrt(422) = rrt(422) * density(25) * density(35) 
  rrt(423) = rrt(423) * density(25) * density(35) 
  rrt(424) = rrt(424) * density(25) * density(35) 
  rrt(425) = rrt(425) * density(24) * density(35) 
  rrt(426) = rrt(426) * density(24) * density(39) 
  rrt(427) = rrt(427) * density(24) * density(39) 
  rrt(428) = rrt(428) * density(24) * density(40) 
  rrt(429) = rrt(429) * density(24) * density(41) 
  rrt(430) = rrt(430) * density(24) * density(37) 
  rrt(431) = rrt(431) * density(25) * density(37) 
  rrt(432) = rrt(432) * density(25) * density(39) 
  rrt(433) = rrt(433) * density(25) * density(40) 
  rrt(434) = rrt(434) * density(25) * density(40) 
  rrt(435) = rrt(435) * density(25) * density(41) 
  rrt(436) = rrt(436) * density(25) * density(41) 
  rrt(437) = rrt(437) * density(25) * density(41) 
  rrt(438) = rrt(438) * density(25) * density(42) 
  rrt(439) = rrt(439) * density(25) * density(42) 
  rrt(440) = rrt(440) * density(25) * density(42) 
  rrt(441) = rrt(441) * density(37) * density(40) 
  rrt(442) = rrt(442) * density(37) * density(40) 
  rrt(443) = rrt(443) * density(36) * density(40) 
  rrt(444) = rrt(444) * density(36) * density(40) 
  rrt(445) = rrt(445) * density(36) * density(40) 
  rrt(446) = rrt(446) * density(40) 
  rrt(447) = rrt(447) * density(36) * density(41) 
  rrt(448) = rrt(448) * density(36) * density(41) 
  rrt(449) = rrt(449) * density(37) * density(41) 
  rrt(450) = rrt(450) * density(41) 
  rrt(451) = rrt(451) * density(37) * density(39) 
  rrt(452) = rrt(452) * density(37) * density(39) 
  rrt(453) = rrt(453) * density(36) * density(39) 
  rrt(454) = rrt(454) * density(36) * density(39) 
  rrt(455) = rrt(455) * density(36) * density(39) 
  rrt(456) = rrt(456) * density(39) 
  rrt(457) = rrt(457) * density(42) 
  rrt(458) = rrt(458) * density(38) 
  rrt(459) = rrt(459) * density(02) * density(03) 
  rrt(460) = rrt(460) * density(03) * density(09) 
  rrt(461) = rrt(461) * density(03) * density(15) 
  rrt(462) = rrt(462) * density(02) * density(04) 
  rrt(463) = rrt(463) * density(04) * density(09) 
  rrt(464) = rrt(464) * density(04) * density(15) 
  rrt(465) = rrt(465) * density(02) * density(04) 
  rrt(466) = rrt(466) * density(04) * density(09) 
  rrt(467) = rrt(467) * density(04) * density(15) 
  rrt(468) = rrt(468) * density(02) * density(05) 
  rrt(469) = rrt(469) * density(05) * density(09) 
  rrt(470) = rrt(470) * density(05) * density(15) 
  rrt(471) = rrt(471) * density(02) * density(05) 
  rrt(472) = rrt(472) * density(05) * density(09) 
  rrt(473) = rrt(473) * density(05) * density(15) 
  rrt(474) = rrt(474) * density(02) * density(05) 
  rrt(475) = rrt(475) * density(02) * density(05) 
  rrt(476) = rrt(476) * density(05) * density(09) 
  rrt(477) = rrt(477) * density(05) * density(15) 
  rrt(478) = rrt(478) * density(02) * density(06) 
  rrt(479) = rrt(479) * density(06) * density(09) 
  rrt(480) = rrt(480) * density(06) * density(15) 
  rrt(481) = rrt(481) * density(02) * density(06) 
  rrt(482) = rrt(482) * density(06) * density(09) 
  rrt(483) = rrt(483) * density(06) * density(15) 
  rrt(484) = rrt(484) * density(02) * density(10) 
  rrt(485) = rrt(485) * density(09) * density(10) 
  rrt(486) = rrt(486) * density(10) * density(15) 
  rrt(487) = rrt(487) * density(02) * density(16) 
  rrt(488) = rrt(488) * density(09) * density(16) 
  rrt(489) = rrt(489) * density(15) * density(16) 
  rrt(490) = rrt(490) * density(02) * density(17) 
  rrt(491) = rrt(491) * density(09) * density(17) 
  rrt(492) = rrt(492) * density(15) * density(17) 
  rrt(493) = rrt(493) * density(02) * density(18) 
  rrt(494) = rrt(494) * density(09) * density(18) 
  rrt(495) = rrt(495) * density(15) * density(18) 
  rrt(496) = rrt(496) * density(07) 
  rrt(497) = rrt(497) * density(08) 
  rrt(498) = rrt(498) * density(11) 
  rrt(499) = rrt(499) * density(12) 
  rrt(500) = rrt(500) * density(13) 
  rrt(501) = rrt(501) * density(14) 
  rrt(502) = rrt(502) * density(19) 
  rrt(503) = rrt(503) * density(20) 
  ydot(01) = +rrt(001)+rrt(002)+rrt(003)+rrt(004)+rrt(005)-rrt(006)+rrt(020)+rrt(021)+rrt(022)+rrt(023)+rrt(024)+rrt(025)+rrt(026)&
             +rrt(027)+rrt(028)+rrt(029)+rrt(030)+rrt(031)+rrt(032)+rrt(033)+rrt(034)+rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)&
             +rrt(048)+rrt(049)+rrt(050)+rrt(051)+rrt(058)+rrt(059)+rrt(060)+rrt(061)+rrt(062)+rrt(063)+rrt(064)+rrt(065)+rrt(066)&
             +rrt(067)+rrt(068)+rrt(069)+rrt(070)-rrt(071)+rrt(083)+rrt(084)+rrt(085)-rrt(086)+rrt(096)+rrt(097)+rrt(098)+rrt(099)&
             +rrt(100)+rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)+rrt(106)+rrt(107)+rrt(108)+rrt(110)+rrt(112)+rrt(113)-rrt(114)&
             +rrt(130)+rrt(131)+rrt(132)+rrt(133)+rrt(134)+rrt(135)-rrt(139)-rrt(140)-rrt(141)+rrt(147)+rrt(148)+rrt(149)+rrt(150)&
             +rrt(152)-rrt(154)-rrt(155)+rrt(156)-rrt(157)-rrt(158)-rrt(159)-rrt(160)-rrt(161)-rrt(162)-rrt(163)-rrt(164)-rrt(165)&
             -rrt(166)-rrt(167)-rrt(168)-rrt(169)-rrt(170)-rrt(171)+rrt(283)+rrt(284)+rrt(291)+rrt(292)+rrt(298)+rrt(299)+rrt(300)&
             +rrt(301)+rrt(302)+rrt(303)+rrt(304)+rrt(305)+rrt(334)+rrt(336)+rrt(337)+rrt(362)+rrt(363)+rrt(364)+rrt(368)+rrt(369)&
             +rrt(370)+rrt(371)+rrt(372)+rrt(373)+rrt(390)+rrt(391)+rrt(394)+rrt(395)+rrt(396)+rrt(397)+rrt(427)+rrt(429)+rrt(432)&
             +rrt(434)+rrt(436)+rrt(446)+rrt(450)+rrt(456) 
  ydot(02) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)-rrt(010)-rrt(011)-rrt(012)-rrt(013)&
             +rrt(014)+rrt(015)+rrt(052)+rrt(053)+rrt(054)+rrt(055)+rrt(159)+rrt(162)+  2.d0 * rrt(163)-rrt(172)-rrt(173)+rrt(182)&
             +rrt(183)+rrt(185)+rrt(186)+rrt(187)+rrt(188)+rrt(189)+rrt(190)+rrt(191)+rrt(192)+rrt(193)+rrt(194)+rrt(197)-rrt(210)&
             -rrt(211)-rrt(212)-rrt(213)-rrt(214)-rrt(215)-rrt(216)-rrt(217)-rrt(218)-rrt(219)+rrt(283)+  2.d0 * rrt(284)+rrt(285)&
             +rrt(286)+rrt(287)+rrt(288)+rrt(291)+  2.d0 * rrt(292)+rrt(293)+rrt(294)+rrt(295)+rrt(296)+rrt(298)+rrt(299)+rrt(300)&
             +rrt(301)+  2.d0 * rrt(302)+  2.d0 * rrt(303)+  2.d0 * rrt(304)+  2.d0 * rrt(305)+rrt(306)+rrt(307)+rrt(308)+rrt(309)&
             +rrt(310)+rrt(311)+rrt(312)+rrt(313)+rrt(314)+rrt(315)+rrt(316)+rrt(317)+rrt(318)+rrt(319)+rrt(320)+rrt(321)+rrt(326)&
             +rrt(338)+rrt(339)+rrt(340)+rrt(374)+rrt(375)+rrt(403)+  2.d0 * rrt(404)+  2.d0 * rrt(405)+rrt(409)+rrt(410)&
             +  2.d0 * rrt(412)+  2.d0 * rrt(413)+rrt(414)+rrt(415)+  3.d0 * rrt(416)+  3.d0 * rrt(417)+  2.d0 * rrt(418)+rrt(419)&
             +rrt(420)+rrt(421)+rrt(423)+rrt(424)+rrt(425)+rrt(459)+rrt(460)+rrt(461)+rrt(462)+rrt(463)+rrt(464)-rrt(474)+rrt(496)&
             +rrt(497) 
  ydot(03) = +rrt(008)-rrt(016)-rrt(020)-rrt(024)-rrt(028)-rrt(032)-rrt(036)-rrt(040)-rrt(044)-rrt(048)-rrt(052)-rrt(174)-rrt(178)&
             -rrt(220)-rrt(224)-rrt(228)-rrt(232)-rrt(236)-rrt(240)-rrt(244)-rrt(248)-rrt(252)-rrt(256)-rrt(459)-rrt(460)-rrt(461)&
             +rrt(465)+rrt(466)+rrt(467)+rrt(474)+rrt(475)+rrt(476)+rrt(477)+rrt(481)+rrt(482)+rrt(483) 
  ydot(04) = +rrt(009)-rrt(017)-rrt(021)-rrt(025)-rrt(029)-rrt(033)-rrt(037)-rrt(041)-rrt(045)-rrt(049)-rrt(053)-rrt(175)-rrt(179)&
             -rrt(221)-rrt(225)-rrt(229)-rrt(233)-rrt(237)-rrt(241)-rrt(245)-rrt(249)-rrt(253)-rrt(257)-rrt(462)-rrt(463)-rrt(464)&
             -rrt(465)-rrt(466)-rrt(467)+rrt(468)+rrt(469)+rrt(470)+rrt(474)+rrt(478)+rrt(479)+rrt(480) 
  ydot(05) = +rrt(010)-rrt(019)-rrt(023)-rrt(027)-rrt(031)-rrt(035)-rrt(039)-rrt(043)-rrt(047)-rrt(051)-rrt(055)-rrt(176)-rrt(180)&
             -rrt(222)-rrt(226)-rrt(230)-rrt(234)-rrt(238)-rrt(242)-rrt(246)-rrt(250)-rrt(254)-rrt(258)-rrt(468)-rrt(469)-rrt(470)&
             -rrt(471)-rrt(472)-rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(477) 
  ydot(06) = +rrt(011)-rrt(018)-rrt(022)-rrt(026)-rrt(030)-rrt(034)-rrt(038)-rrt(042)-rrt(046)-rrt(050)-rrt(054)-rrt(177)-rrt(181)&
             -rrt(223)-rrt(227)-rrt(231)-rrt(235)-rrt(239)-rrt(243)-rrt(247)-rrt(251)-rrt(255)-rrt(259)+rrt(471)+rrt(472)+rrt(473)&
             -rrt(478)-rrt(479)-rrt(480)-rrt(481)-rrt(482)-rrt(483) 
  ydot(07) = +rrt(012)-rrt(014)+rrt(040)+rrt(041)+rrt(042)+rrt(043)-rrt(056)-rrt(058)-rrt(060)-rrt(062)-rrt(064)-rrt(066)-rrt(260)&
             -rrt(262)-rrt(264)-rrt(266)-rrt(268)-rrt(270)-rrt(272)-rrt(274)-rrt(276)-rrt(278)-rrt(280)-rrt(496) 
  ydot(08) = +rrt(013)-rrt(015)+rrt(044)+rrt(045)+rrt(046)+rrt(047)-rrt(057)-rrt(059)-rrt(061)-rrt(063)-rrt(065)-rrt(067)-rrt(261)&
             -rrt(263)-rrt(265)-rrt(267)-rrt(269)-rrt(271)-rrt(273)-rrt(275)-rrt(277)-rrt(279)-rrt(281)-rrt(497) 
  ydot(09) = +rrt(004)+rrt(006)+rrt(007)+rrt(016)+rrt(017)+rrt(018)+rrt(019)+rrt(032)+rrt(033)+rrt(034)+rrt(035)+rrt(036)+rrt(037)&
             +rrt(038)+rrt(039)+rrt(056)+rrt(057)+rrt(064)+rrt(065)-rrt(068)-rrt(069)-rrt(070)-rrt(071)-rrt(072)-rrt(073)-rrt(074)&
             -rrt(075)-rrt(076)-rrt(077)+rrt(078)+rrt(079)+rrt(080)+rrt(081)+rrt(091)+rrt(157)+  2.d0 * rrt(161)+rrt(162)+rrt(172)&
             +  2.d0 * rrt(173)+rrt(174)+rrt(175)+rrt(176)+rrt(177)+  2.d0 * rrt(178)+  2.d0 * rrt(179)+  2.d0 * rrt(180)&
             +  2.d0 * rrt(181)-rrt(182)-rrt(183)-rrt(184)-rrt(185)-rrt(187)-rrt(188)-rrt(189)+rrt(196)+rrt(197)+rrt(199)+rrt(200)&
             +rrt(201)+rrt(205)+  2.d0 * rrt(206)+rrt(211)+rrt(213)+rrt(214)+rrt(220)+rrt(221)+rrt(222)+rrt(223)+rrt(228)+rrt(229)&
             +rrt(230)+rrt(231)+rrt(232)+rrt(233)+rrt(234)+rrt(235)+rrt(260)+rrt(261)+rrt(264)+rrt(265)+rrt(266)+rrt(267)-rrt(282)&
             -rrt(283)-rrt(284)-rrt(285)-rrt(286)-rrt(287)-rrt(288)-rrt(289)+rrt(327)+  2.d0 * rrt(328)+rrt(329)+rrt(341)+rrt(342)&
             +rrt(343)+  2.d0 * rrt(344)+  2.d0 * rrt(345)+  2.d0 * rrt(346)+rrt(347)+rrt(348)+rrt(349)+rrt(376)+rrt(377)&
             +  2.d0 * rrt(378)+  2.d0 * rrt(379)+rrt(380)+rrt(381)+rrt(398)+rrt(402)+rrt(406)+rrt(407)+rrt(408)+  2.d0 * rrt(409)&
             +  2.d0 * rrt(410)+  2.d0 * rrt(411)+rrt(412)+rrt(413)+rrt(414)+rrt(484)+rrt(485)+rrt(486)+rrt(498)+rrt(499)+rrt(500)&
             +rrt(501) 
  ydot(10) = +rrt(073)-rrt(082)-rrt(083)-rrt(084)-rrt(085)-rrt(086)-rrt(087)-rrt(088)-rrt(089)-rrt(090)-rrt(091)-rrt(186)-rrt(190)&
             -rrt(191)-rrt(192)-rrt(193)-rrt(194)-rrt(195)-rrt(290)-rrt(291)-rrt(292)-rrt(293)-rrt(294)-rrt(295)-rrt(296)-rrt(297)&
             -rrt(484)-rrt(485)-rrt(486) 
  ydot(11) = +rrt(074)-rrt(078)+rrt(087)-rrt(092)-rrt(096)-rrt(100)-rrt(104)-rrt(298)-rrt(302)-rrt(306)-rrt(310)-rrt(314)-rrt(318)&
             -rrt(322)-rrt(498) 
  ydot(12) = +rrt(075)-rrt(079)+rrt(088)-rrt(093)-rrt(097)-rrt(101)-rrt(105)-rrt(299)-rrt(303)-rrt(307)-rrt(311)-rrt(315)-rrt(319)&
             -rrt(323)-rrt(499) 
  ydot(13) = +rrt(076)-rrt(080)+rrt(089)-rrt(094)-rrt(098)-rrt(102)-rrt(106)-rrt(300)-rrt(304)-rrt(308)-rrt(312)-rrt(316)-rrt(320)&
             -rrt(324)-rrt(500) 
  ydot(14) = +rrt(077)-rrt(081)+rrt(090)-rrt(095)-rrt(099)-rrt(103)-rrt(107)-rrt(301)-rrt(305)-rrt(309)-rrt(313)-rrt(317)-rrt(321)&
             -rrt(325)-rrt(501) 
  ydot(15) = +rrt(003)+rrt(028)+rrt(029)+rrt(030)+rrt(031)+rrt(062)+rrt(063)-rrt(111)-rrt(112)-rrt(113)-rrt(114)-rrt(115)-rrt(116)&
             -rrt(117)-rrt(118)-rrt(119)+rrt(120)+rrt(121)+rrt(122)+rrt(123)+rrt(142)+rrt(143)+rrt(144)+rrt(151)+rrt(154)+rrt(158)&
             +rrt(159)-rrt(165)+rrt(168)+  2.d0 * rrt(171)+rrt(172)+rrt(174)+rrt(175)+rrt(176)+rrt(177)+rrt(183)-rrt(185)-rrt(190)&
             +rrt(194)-rrt(196)-rrt(197)-rrt(198)+  2.d0 * rrt(207)+rrt(208)+rrt(209)+rrt(217)+rrt(218)+rrt(244)+rrt(245)+rrt(246)&
             +rrt(247)+rrt(248)+rrt(249)+rrt(250)+rrt(251)+rrt(272)+rrt(273)+rrt(274)+rrt(275)-rrt(326)-rrt(327)-rrt(328)-rrt(329)&
             -rrt(330)-rrt(331)-rrt(332)-rrt(333)-rrt(334)-rrt(335)+rrt(336)+rrt(368)+rrt(369)+rrt(370)+rrt(394)+rrt(395)+rrt(400)&
             +rrt(405)+rrt(406)+rrt(410)+rrt(411)+rrt(413)+rrt(414)+rrt(417)+rrt(418)+rrt(419)+  2.d0 * rrt(421)+rrt(422)+rrt(423)&
             +rrt(425)+  2.d0 * rrt(427)+rrt(428)+  3.d0 * rrt(429)+rrt(430)+rrt(432)+rrt(433)+  2.d0 * rrt(436)+rrt(437)+rrt(438)&
             +  2.d0 * rrt(439)+rrt(442)+  2.d0 * rrt(443)+rrt(444)+  2.d0 * rrt(445)+rrt(446)+rrt(447)+rrt(452)+rrt(453)+rrt(457)&
             +rrt(458)+rrt(487)+rrt(488)+rrt(489)+rrt(490)+rrt(491)+rrt(492)+rrt(493)+rrt(494)+rrt(495)+rrt(502)+rrt(503) 
  ydot(16) = +rrt(115)-rrt(124)-rrt(127)-rrt(130)-rrt(133)-rrt(136)-rrt(139)-rrt(142)-rrt(187)-rrt(191)-rrt(199)-rrt(202)-rrt(338)&
             -rrt(341)-rrt(344)-rrt(347)-rrt(350)-rrt(353)-rrt(356)-rrt(359)-rrt(362)-rrt(365)-rrt(487)-rrt(488)-rrt(489) 
  ydot(17) = +rrt(116)-rrt(125)-rrt(128)-rrt(131)-rrt(134)-rrt(137)-rrt(140)-rrt(143)-rrt(188)-rrt(192)-rrt(200)-rrt(203)-rrt(339)&
             -rrt(342)-rrt(345)-rrt(348)-rrt(351)-rrt(354)-rrt(357)-rrt(360)-rrt(363)-rrt(366)-rrt(490)-rrt(491)-rrt(492) 
  ydot(18) = +rrt(117)-rrt(126)-rrt(129)-rrt(132)-rrt(135)-rrt(138)-rrt(141)-rrt(144)-rrt(189)-rrt(193)-rrt(201)-rrt(204)-rrt(340)&
             -rrt(343)-rrt(346)-rrt(349)-rrt(352)-rrt(355)-rrt(358)-rrt(361)-rrt(364)-rrt(367)-rrt(493)-rrt(494)-rrt(495) 
  ydot(19) = +rrt(118)-rrt(120)-rrt(122)+rrt(124)+rrt(125)+rrt(126)-rrt(145)-rrt(147)-rrt(149)-rrt(374)-rrt(376)-rrt(378)-rrt(380)&
             -rrt(382)-rrt(384)-rrt(386)-rrt(388)-rrt(390)-rrt(392)-rrt(502) 
  ydot(20) = +rrt(119)-rrt(121)-rrt(123)+rrt(127)+rrt(128)+rrt(129)-rrt(146)-rrt(148)-rrt(150)-rrt(375)-rrt(377)-rrt(379)-rrt(381)&
             -rrt(383)-rrt(385)-rrt(387)-rrt(389)-rrt(391)-rrt(393)-rrt(503) 
  ydot(21) = +rrt(184)+rrt(195)-rrt(197)-rrt(206) 
  ydot(22) = +rrt(005)+rrt(048)+rrt(049)+rrt(050)+rrt(051)+rrt(066)+rrt(067)+rrt(070)+rrt(071)+rrt(072)+rrt(082)+rrt(085)+rrt(086)&
             +rrt(092)+rrt(093)+rrt(094)+rrt(095)+rrt(104)+rrt(105)+rrt(106)+rrt(107)-rrt(108)+  2.d0 * rrt(109)+rrt(158)+rrt(160)&
             +  2.d0 * rrt(164)-rrt(173)-rrt(178)-rrt(179)-rrt(180)-rrt(181)-rrt(184)-rrt(195)-rrt(196)-rrt(199)-rrt(200)-rrt(201)&
             -rrt(205)+rrt(289)+rrt(297)+rrt(322)+rrt(323)+rrt(324)+rrt(325)-rrt(398)-rrt(399)-rrt(400)-rrt(401) 
  ydot(23) = -rrt(109)-rrt(110)+rrt(401) 
  ydot(24) = -rrt(151)-rrt(152)-rrt(153)-rrt(154)-rrt(155)-rrt(166)-rrt(183)-rrt(194)+rrt(198)+rrt(202)+rrt(203)+rrt(204)-rrt(207)&
             -rrt(208)+rrt(334)+rrt(337)+rrt(362)+rrt(363)+rrt(364)+rrt(371)+rrt(372)+rrt(373)+rrt(390)+rrt(391)+rrt(396)+rrt(397)&
             -rrt(425)-rrt(426)-rrt(427)-rrt(428)-rrt(429)-rrt(430)+rrt(434)+rrt(435)+rrt(440)+rrt(441)+rrt(447)+rrt(448)+rrt(449)&
             +rrt(450)+rrt(455) 
  ydot(25) = +rrt(002)+rrt(007)+rrt(016)+rrt(017)+rrt(018)+rrt(019)+rrt(024)+rrt(025)+rrt(026)+rrt(027)+rrt(056)+rrt(057)+rrt(060)&
             +rrt(061)+rrt(069)+rrt(072)+rrt(082)+rrt(084)+rrt(092)+rrt(093)+rrt(094)+rrt(095)+rrt(100)+rrt(101)+rrt(102)+rrt(103)&
             +  2.d0 * rrt(111)+rrt(113)+rrt(114)+rrt(133)+rrt(134)+rrt(135)+  2.d0 * rrt(136)+  2.d0 * rrt(137)+  2.d0 * rrt(138)&
             +rrt(139)+rrt(140)+rrt(141)+  2.d0 * rrt(145)+  2.d0 * rrt(146)+rrt(149)+rrt(150)+rrt(151)+rrt(152)+rrt(153)+rrt(155)&
             -rrt(156)+rrt(157)+rrt(160)-rrt(167)+  2.d0 * rrt(169)+rrt(170)-rrt(172)-rrt(174)-rrt(175)-rrt(176)-rrt(177)-rrt(182)&
             +rrt(185)-rrt(186)+rrt(187)+rrt(188)+rrt(189)+rrt(190)+rrt(191)+rrt(192)+rrt(193)+rrt(196)-rrt(198)+rrt(199)+rrt(200)&
             +rrt(201)-rrt(202)-rrt(203)-rrt(204)-rrt(205)-rrt(206)-rrt(207)+rrt(208)-  2.d0 * rrt(209)+rrt(212)+rrt(224)+rrt(225)&
             +rrt(226)+rrt(227)+rrt(262)+rrt(263)+rrt(282)+rrt(290)+rrt(330)+rrt(331)+rrt(350)+rrt(351)+rrt(352)+rrt(353)+rrt(354)&
             +rrt(355)+rrt(382)+rrt(383)+rrt(384)+rrt(385)+rrt(399)-rrt(402)-rrt(403)+rrt(404)+rrt(406)-rrt(407)+rrt(409)+rrt(412)&
             +rrt(416)+rrt(419)-rrt(420)-rrt(422)-rrt(423)-rrt(424)+rrt(426)-rrt(431)-rrt(432)-rrt(433)-rrt(434)-rrt(435)-rrt(436)&
             -rrt(437)-rrt(438)-rrt(439)-rrt(440)+rrt(442)+  2.d0 * rrt(444)+  2.d0 * rrt(448)+rrt(449)+  2.d0 * rrt(451)+rrt(453)&
             +  3.d0 * rrt(454)+rrt(456) 
  ydot(26) = +rrt(001)+rrt(020)+rrt(021)+rrt(022)+rrt(023)+rrt(058)+rrt(059)-rrt(157)-rrt(158)+rrt(212)+rrt(214)-rrt(219)+rrt(224)&
             +rrt(225)+rrt(226)+rrt(227)+rrt(232)+rrt(233)+rrt(234)+rrt(235)-rrt(252)-rrt(253)-rrt(254)-rrt(255)+rrt(262)+rrt(263)&
             +rrt(266)+rrt(267)-rrt(276)-rrt(277)-rrt(326)-rrt(338)-rrt(339)-rrt(340)-rrt(374)-rrt(375)-rrt(402)-rrt(403)-rrt(404)&
             -rrt(405)-rrt(406)+rrt(415) 
  ydot(27) = -rrt(159)+rrt(210)+rrt(256)+rrt(257)+rrt(258)+rrt(259)+rrt(278)+rrt(279)+rrt(280)+rrt(281) 
  ydot(28) = +rrt(002)+rrt(024)+rrt(025)+rrt(026)+rrt(027)+rrt(060)+rrt(061)+rrt(068)+rrt(083)+rrt(096)+rrt(097)+rrt(098)+rrt(099)&
             -rrt(160)+rrt(213)-rrt(214)+rrt(228)+rrt(229)+rrt(230)+rrt(231)-rrt(232)-rrt(233)-rrt(234)-rrt(235)+rrt(264)+rrt(265)&
             -rrt(266)-rrt(267)+rrt(282)+rrt(289)+rrt(290)+rrt(297)+rrt(322)+rrt(323)+rrt(324)+rrt(325)-rrt(327)+rrt(330)-rrt(341)&
             -rrt(342)-rrt(343)+rrt(350)+rrt(351)+rrt(352)-rrt(376)-rrt(377)+rrt(382)+rrt(383)-rrt(398)+rrt(399)-rrt(407)+rrt(408) 
  ydot(29) = -rrt(161)+rrt(285)+rrt(287)+rrt(293)+rrt(295)+rrt(306)+rrt(307)+rrt(308)+rrt(309)+rrt(314)+rrt(315)+rrt(316)+rrt(317)&
             -rrt(328)-rrt(344)-rrt(345)-rrt(346)-rrt(378)-rrt(379)-rrt(408)-rrt(409)-rrt(410)-rrt(411) 
  ydot(30) = -rrt(162)-rrt(285)+rrt(286)-rrt(287)+rrt(288)-rrt(293)+rrt(294)-rrt(295)+rrt(296)-rrt(306)-rrt(307)-rrt(308)-rrt(309)&
             +rrt(310)+rrt(311)+rrt(312)+rrt(313)-rrt(314)-rrt(315)-rrt(316)-rrt(317)+rrt(318)+rrt(319)+rrt(320)+rrt(321)-rrt(412)&
             -rrt(413)-rrt(414) 
  ydot(31) = -rrt(163)+rrt(219)+rrt(252)+rrt(253)+rrt(254)+rrt(255)+rrt(276)+rrt(277)-rrt(286)-rrt(288)-rrt(294)-rrt(296)-rrt(310)&
             -rrt(311)-rrt(312)-rrt(313)-rrt(318)-rrt(319)-rrt(320)-rrt(321)-rrt(415)-rrt(416)-rrt(417)-rrt(418) 
  ydot(32) = +rrt(003)+rrt(028)+rrt(029)+rrt(030)+rrt(031)+rrt(062)+rrt(063)+rrt(069)+rrt(084)+rrt(100)+rrt(101)+rrt(102)+rrt(103)&
             +rrt(108)-rrt(213)-rrt(228)-rrt(229)-rrt(230)-rrt(231)-rrt(264)-rrt(265)-rrt(289)-rrt(297)-rrt(322)-rrt(323)-rrt(324)&
             -rrt(325)-rrt(329)-rrt(330)-rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(380)-rrt(381)-rrt(382)-rrt(383)&
             +rrt(398)+rrt(400)+rrt(401) 
  ydot(33) = +rrt(110)-rrt(164)-rrt(401) 
  ydot(34) = +rrt(215)+rrt(217)+rrt(236)+rrt(237)+rrt(238)+rrt(239)+rrt(244)+rrt(245)+rrt(246)+rrt(247)+rrt(268)+rrt(269)+rrt(272)&
             +rrt(273)-rrt(284)-rrt(292)-rrt(302)-rrt(303)-rrt(304)-rrt(305)-rrt(404)-rrt(409)-rrt(412)-rrt(416)-rrt(419)-rrt(420)&
             +rrt(422) 
  ydot(35) = +rrt(216)+rrt(218)+rrt(240)+rrt(241)+rrt(242)+rrt(243)+rrt(248)+rrt(249)+rrt(250)+rrt(251)+rrt(270)+rrt(271)+rrt(274)&
             +rrt(275)-rrt(405)-rrt(410)-rrt(413)-rrt(417)-rrt(421)-rrt(422)-rrt(423)-rrt(424)-rrt(425) 
  ydot(36) = +rrt(005)+rrt(048)+rrt(049)+rrt(050)+rrt(051)+rrt(066)+rrt(067)+rrt(112)+rrt(130)+rrt(131)+rrt(132)+rrt(147)+rrt(148)&
             +rrt(152)-rrt(168)-rrt(169)-rrt(210)+rrt(211)+rrt(220)+rrt(221)+rrt(222)+rrt(223)-rrt(256)-rrt(257)-rrt(258)-rrt(259)&
             +rrt(260)+rrt(261)-rrt(278)-rrt(279)-rrt(280)-rrt(281)+rrt(326)+rrt(327)+rrt(328)+rrt(331)-rrt(332)+rrt(338)+rrt(339)&
             +rrt(340)+rrt(341)+rrt(342)+rrt(343)+rrt(344)+rrt(345)+rrt(346)+rrt(353)+rrt(354)+rrt(355)-rrt(356)-rrt(357)-rrt(358)&
             +rrt(374)+rrt(375)+rrt(376)+rrt(377)+rrt(378)+rrt(379)+rrt(384)+rrt(385)-rrt(386)-rrt(387)-rrt(399)-rrt(400)+rrt(402)&
             -rrt(419)-rrt(421)+rrt(430)+rrt(431)-rrt(443)-rrt(444)-rrt(445)-rrt(447)-rrt(448)-rrt(453)-rrt(454)-rrt(455)+rrt(458) 
  ydot(37) = +rrt(004)+rrt(032)+rrt(033)+rrt(034)+rrt(035)+rrt(064)+rrt(065)+rrt(070)+rrt(085)+rrt(104)+rrt(105)+rrt(106)+rrt(107)&
             +rrt(113)+rrt(133)+rrt(134)+rrt(135)+rrt(149)+rrt(150)+rrt(153)+rrt(156)-rrt(170)-rrt(211)-rrt(212)-rrt(220)-rrt(221)&
             -rrt(222)-rrt(223)-rrt(224)-rrt(225)-rrt(226)-rrt(227)-rrt(260)-rrt(261)-rrt(262)-rrt(263)-rrt(282)-rrt(290)+rrt(329)&
             -rrt(331)+rrt(347)+rrt(348)+rrt(349)-rrt(353)-rrt(354)-rrt(355)+rrt(380)+rrt(381)-rrt(384)-rrt(385)+rrt(403)+rrt(407)&
             -rrt(430)-rrt(431)-rrt(441)-rrt(442)-rrt(449)-rrt(451)-rrt(452) 
  ydot(38) = -rrt(171)+rrt(332)+rrt(356)+rrt(357)+rrt(358)+rrt(386)+rrt(387)-rrt(458) 
  ydot(39) = +rrt(006)+rrt(036)+rrt(037)+rrt(038)+rrt(039)+rrt(071)+rrt(086)+rrt(114)+rrt(139)+rrt(140)+rrt(141)+rrt(153)+rrt(154)&
             +rrt(167)-rrt(215)-rrt(236)-rrt(237)-rrt(238)-rrt(239)-rrt(268)-rrt(269)-rrt(283)-rrt(291)-rrt(298)-rrt(299)-rrt(300)&
             -rrt(301)-rrt(334)-rrt(335)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)-rrt(367)-rrt(390)-rrt(391)-rrt(392)-rrt(393)&
             +rrt(423)-rrt(426)-rrt(427)-rrt(432)+rrt(433)+rrt(435)+rrt(439)-rrt(451)-rrt(452)-rrt(453)-rrt(454)-rrt(455)-rrt(456) 
  ydot(40) = +rrt(155)+rrt(165)-rrt(216)-rrt(240)-rrt(241)-rrt(242)-rrt(243)-rrt(270)-rrt(271)-rrt(333)-rrt(336)-rrt(359)-rrt(360)&
             -rrt(361)-rrt(368)-rrt(369)-rrt(370)-rrt(388)-rrt(389)-rrt(394)-rrt(395)-rrt(406)-rrt(411)-rrt(414)-rrt(418)+rrt(420)&
             -rrt(428)-rrt(433)-rrt(434)+rrt(437)+rrt(440)-rrt(441)-rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)+rrt(457) 
  ydot(41) = +rrt(166)-rrt(217)-rrt(244)-rrt(245)-rrt(246)-rrt(247)-rrt(272)-rrt(273)+rrt(335)-rrt(337)+rrt(365)+rrt(366)+rrt(367)&
             -rrt(371)-rrt(372)-rrt(373)+rrt(392)+rrt(393)-rrt(396)-rrt(397)+rrt(424)+rrt(425)+rrt(426)+rrt(428)-rrt(429)-rrt(435)&
             -rrt(436)-rrt(437)+rrt(438)-rrt(447)-rrt(448)-rrt(449)-rrt(450) 
  ydot(42) = -rrt(218)-rrt(248)-rrt(249)-rrt(250)-rrt(251)-rrt(274)-rrt(275)+rrt(333)+rrt(359)+rrt(360)+rrt(361)+rrt(388)+rrt(389)&
             -rrt(438)-rrt(439)-rrt(440)-rrt(457) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(43) = 0.0d0
  if( lgas_heating ) then
    ydot(43) = ( ZDPlasKin_cfg(14)/k_B + ydot(43) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(43) = ydot(43) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(43)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) + rrt(001) * density(02) 
  pd(01,02) = pd(01,02) + rrt(001) * density(01) 
  pd(02,01) = pd(02,01) - rrt(001) * density(02) 
  pd(02,02) = pd(02,02) - rrt(001) * density(01) 
  pd(26,01) = pd(26,01) + rrt(001) * density(02) 
  pd(26,02) = pd(26,02) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) + rrt(002) * density(02) 
  pd(01,02) = pd(01,02) + rrt(002) * density(01) 
  pd(02,01) = pd(02,01) - rrt(002) * density(02) 
  pd(02,02) = pd(02,02) - rrt(002) * density(01) 
  pd(25,01) = pd(25,01) + rrt(002) * density(02) 
  pd(25,02) = pd(25,02) + rrt(002) * density(01) 
  pd(28,01) = pd(28,01) + rrt(002) * density(02) 
  pd(28,02) = pd(28,02) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) + rrt(003) * density(02) 
  pd(01,02) = pd(01,02) + rrt(003) * density(01) 
  pd(02,01) = pd(02,01) - rrt(003) * density(02) 
  pd(02,02) = pd(02,02) - rrt(003) * density(01) 
  pd(15,01) = pd(15,01) + rrt(003) * density(02) 
  pd(15,02) = pd(15,02) + rrt(003) * density(01) 
  pd(32,01) = pd(32,01) + rrt(003) * density(02) 
  pd(32,02) = pd(32,02) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) + rrt(004) * density(02) 
  pd(01,02) = pd(01,02) + rrt(004) * density(01) 
  pd(02,01) = pd(02,01) - rrt(004) * density(02) 
  pd(02,02) = pd(02,02) - rrt(004) * density(01) 
  pd(09,01) = pd(09,01) + rrt(004) * density(02) 
  pd(09,02) = pd(09,02) + rrt(004) * density(01) 
  pd(37,01) = pd(37,01) + rrt(004) * density(02) 
  pd(37,02) = pd(37,02) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) + rrt(005) * density(02) 
  pd(01,02) = pd(01,02) + rrt(005) * density(01) 
  pd(02,01) = pd(02,01) - rrt(005) * density(02) 
  pd(02,02) = pd(02,02) - rrt(005) * density(01) 
  pd(22,01) = pd(22,01) + rrt(005) * density(02) 
  pd(22,02) = pd(22,02) + rrt(005) * density(01) 
  pd(36,01) = pd(36,01) + rrt(005) * density(02) 
  pd(36,02) = pd(36,02) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(02) 
  pd(01,02) = pd(01,02) - rrt(006) * density(01) 
  pd(02,01) = pd(02,01) - rrt(006) * density(02) 
  pd(02,02) = pd(02,02) - rrt(006) * density(01) 
  pd(09,01) = pd(09,01) + rrt(006) * density(02) 
  pd(09,02) = pd(09,02) + rrt(006) * density(01) 
  pd(39,01) = pd(39,01) + rrt(006) * density(02) 
  pd(39,02) = pd(39,02) + rrt(006) * density(01) 
  pd(02,01) = pd(02,01) - rrt(007) * density(02) 
  pd(02,02) = pd(02,02) - rrt(007) * density(01) 
  pd(09,01) = pd(09,01) + rrt(007) * density(02) 
  pd(09,02) = pd(09,02) + rrt(007) * density(01) 
  pd(25,01) = pd(25,01) + rrt(007) * density(02) 
  pd(25,02) = pd(25,02) + rrt(007) * density(01) 
  pd(02,01) = pd(02,01) - rrt(008) * density(02) 
  pd(02,02) = pd(02,02) - rrt(008) * density(01) 
  pd(03,01) = pd(03,01) + rrt(008) * density(02) 
  pd(03,02) = pd(03,02) + rrt(008) * density(01) 
  pd(02,01) = pd(02,01) - rrt(009) * density(02) 
  pd(02,02) = pd(02,02) - rrt(009) * density(01) 
  pd(04,01) = pd(04,01) + rrt(009) * density(02) 
  pd(04,02) = pd(04,02) + rrt(009) * density(01) 
  pd(02,01) = pd(02,01) - rrt(010) * density(02) 
  pd(02,02) = pd(02,02) - rrt(010) * density(01) 
  pd(05,01) = pd(05,01) + rrt(010) * density(02) 
  pd(05,02) = pd(05,02) + rrt(010) * density(01) 
  pd(02,01) = pd(02,01) - rrt(011) * density(02) 
  pd(02,02) = pd(02,02) - rrt(011) * density(01) 
  pd(06,01) = pd(06,01) + rrt(011) * density(02) 
  pd(06,02) = pd(06,02) + rrt(011) * density(01) 
  pd(02,01) = pd(02,01) - rrt(012) * density(02) 
  pd(02,02) = pd(02,02) - rrt(012) * density(01) 
  pd(07,01) = pd(07,01) + rrt(012) * density(02) 
  pd(07,02) = pd(07,02) + rrt(012) * density(01) 
  pd(02,01) = pd(02,01) - rrt(013) * density(02) 
  pd(02,02) = pd(02,02) - rrt(013) * density(01) 
  pd(08,01) = pd(08,01) + rrt(013) * density(02) 
  pd(08,02) = pd(08,02) + rrt(013) * density(01) 
  pd(02,01) = pd(02,01) + rrt(014) * density(07) 
  pd(02,07) = pd(02,07) + rrt(014) * density(01) 
  pd(07,01) = pd(07,01) - rrt(014) * density(07) 
  pd(07,07) = pd(07,07) - rrt(014) * density(01) 
  pd(02,01) = pd(02,01) + rrt(015) * density(08) 
  pd(02,08) = pd(02,08) + rrt(015) * density(01) 
  pd(08,01) = pd(08,01) - rrt(015) * density(08) 
  pd(08,08) = pd(08,08) - rrt(015) * density(01) 
  pd(03,01) = pd(03,01) - rrt(016) * density(03) 
  pd(03,03) = pd(03,03) - rrt(016) * density(01) 
  pd(09,01) = pd(09,01) + rrt(016) * density(03) 
  pd(09,03) = pd(09,03) + rrt(016) * density(01) 
  pd(25,01) = pd(25,01) + rrt(016) * density(03) 
  pd(25,03) = pd(25,03) + rrt(016) * density(01) 
  pd(04,01) = pd(04,01) - rrt(017) * density(04) 
  pd(04,04) = pd(04,04) - rrt(017) * density(01) 
  pd(09,01) = pd(09,01) + rrt(017) * density(04) 
  pd(09,04) = pd(09,04) + rrt(017) * density(01) 
  pd(25,01) = pd(25,01) + rrt(017) * density(04) 
  pd(25,04) = pd(25,04) + rrt(017) * density(01) 
  pd(06,01) = pd(06,01) - rrt(018) * density(06) 
  pd(06,06) = pd(06,06) - rrt(018) * density(01) 
  pd(09,01) = pd(09,01) + rrt(018) * density(06) 
  pd(09,06) = pd(09,06) + rrt(018) * density(01) 
  pd(25,01) = pd(25,01) + rrt(018) * density(06) 
  pd(25,06) = pd(25,06) + rrt(018) * density(01) 
  pd(05,01) = pd(05,01) - rrt(019) * density(05) 
  pd(05,05) = pd(05,05) - rrt(019) * density(01) 
  pd(09,01) = pd(09,01) + rrt(019) * density(05) 
  pd(09,05) = pd(09,05) + rrt(019) * density(01) 
  pd(25,01) = pd(25,01) + rrt(019) * density(05) 
  pd(25,05) = pd(25,05) + rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(03) 
  pd(01,03) = pd(01,03) + rrt(020) * density(01) 
  pd(03,01) = pd(03,01) - rrt(020) * density(03) 
  pd(03,03) = pd(03,03) - rrt(020) * density(01) 
  pd(26,01) = pd(26,01) + rrt(020) * density(03) 
  pd(26,03) = pd(26,03) + rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(04) 
  pd(01,04) = pd(01,04) + rrt(021) * density(01) 
  pd(04,01) = pd(04,01) - rrt(021) * density(04) 
  pd(04,04) = pd(04,04) - rrt(021) * density(01) 
  pd(26,01) = pd(26,01) + rrt(021) * density(04) 
  pd(26,04) = pd(26,04) + rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(06) 
  pd(01,06) = pd(01,06) + rrt(022) * density(01) 
  pd(06,01) = pd(06,01) - rrt(022) * density(06) 
  pd(06,06) = pd(06,06) - rrt(022) * density(01) 
  pd(26,01) = pd(26,01) + rrt(022) * density(06) 
  pd(26,06) = pd(26,06) + rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(05) 
  pd(01,05) = pd(01,05) + rrt(023) * density(01) 
  pd(05,01) = pd(05,01) - rrt(023) * density(05) 
  pd(05,05) = pd(05,05) - rrt(023) * density(01) 
  pd(26,01) = pd(26,01) + rrt(023) * density(05) 
  pd(26,05) = pd(26,05) + rrt(023) * density(01) 
  pd(01,01) = pd(01,01) + rrt(024) * density(03) 
  pd(01,03) = pd(01,03) + rrt(024) * density(01) 
  pd(03,01) = pd(03,01) - rrt(024) * density(03) 
  pd(03,03) = pd(03,03) - rrt(024) * density(01) 
  pd(25,01) = pd(25,01) + rrt(024) * density(03) 
  pd(25,03) = pd(25,03) + rrt(024) * density(01) 
  pd(28,01) = pd(28,01) + rrt(024) * density(03) 
  pd(28,03) = pd(28,03) + rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(04) 
  pd(01,04) = pd(01,04) + rrt(025) * density(01) 
  pd(04,01) = pd(04,01) - rrt(025) * density(04) 
  pd(04,04) = pd(04,04) - rrt(025) * density(01) 
  pd(25,01) = pd(25,01) + rrt(025) * density(04) 
  pd(25,04) = pd(25,04) + rrt(025) * density(01) 
  pd(28,01) = pd(28,01) + rrt(025) * density(04) 
  pd(28,04) = pd(28,04) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(06) 
  pd(01,06) = pd(01,06) + rrt(026) * density(01) 
  pd(06,01) = pd(06,01) - rrt(026) * density(06) 
  pd(06,06) = pd(06,06) - rrt(026) * density(01) 
  pd(25,01) = pd(25,01) + rrt(026) * density(06) 
  pd(25,06) = pd(25,06) + rrt(026) * density(01) 
  pd(28,01) = pd(28,01) + rrt(026) * density(06) 
  pd(28,06) = pd(28,06) + rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(05) 
  pd(01,05) = pd(01,05) + rrt(027) * density(01) 
  pd(05,01) = pd(05,01) - rrt(027) * density(05) 
  pd(05,05) = pd(05,05) - rrt(027) * density(01) 
  pd(25,01) = pd(25,01) + rrt(027) * density(05) 
  pd(25,05) = pd(25,05) + rrt(027) * density(01) 
  pd(28,01) = pd(28,01) + rrt(027) * density(05) 
  pd(28,05) = pd(28,05) + rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(03) 
  pd(01,03) = pd(01,03) + rrt(028) * density(01) 
  pd(03,01) = pd(03,01) - rrt(028) * density(03) 
  pd(03,03) = pd(03,03) - rrt(028) * density(01) 
  pd(15,01) = pd(15,01) + rrt(028) * density(03) 
  pd(15,03) = pd(15,03) + rrt(028) * density(01) 
  pd(32,01) = pd(32,01) + rrt(028) * density(03) 
  pd(32,03) = pd(32,03) + rrt(028) * density(01) 
  pd(01,01) = pd(01,01) + rrt(029) * density(04) 
  pd(01,04) = pd(01,04) + rrt(029) * density(01) 
  pd(04,01) = pd(04,01) - rrt(029) * density(04) 
  pd(04,04) = pd(04,04) - rrt(029) * density(01) 
  pd(15,01) = pd(15,01) + rrt(029) * density(04) 
  pd(15,04) = pd(15,04) + rrt(029) * density(01) 
  pd(32,01) = pd(32,01) + rrt(029) * density(04) 
  pd(32,04) = pd(32,04) + rrt(029) * density(01) 
  pd(01,01) = pd(01,01) + rrt(030) * density(06) 
  pd(01,06) = pd(01,06) + rrt(030) * density(01) 
  pd(06,01) = pd(06,01) - rrt(030) * density(06) 
  pd(06,06) = pd(06,06) - rrt(030) * density(01) 
  pd(15,01) = pd(15,01) + rrt(030) * density(06) 
  pd(15,06) = pd(15,06) + rrt(030) * density(01) 
  pd(32,01) = pd(32,01) + rrt(030) * density(06) 
  pd(32,06) = pd(32,06) + rrt(030) * density(01) 
  pd(01,01) = pd(01,01) + rrt(031) * density(05) 
  pd(01,05) = pd(01,05) + rrt(031) * density(01) 
  pd(05,01) = pd(05,01) - rrt(031) * density(05) 
  pd(05,05) = pd(05,05) - rrt(031) * density(01) 
  pd(15,01) = pd(15,01) + rrt(031) * density(05) 
  pd(15,05) = pd(15,05) + rrt(031) * density(01) 
  pd(32,01) = pd(32,01) + rrt(031) * density(05) 
  pd(32,05) = pd(32,05) + rrt(031) * density(01) 
  pd(01,01) = pd(01,01) + rrt(032) * density(03) 
  pd(01,03) = pd(01,03) + rrt(032) * density(01) 
  pd(03,01) = pd(03,01) - rrt(032) * density(03) 
  pd(03,03) = pd(03,03) - rrt(032) * density(01) 
  pd(09,01) = pd(09,01) + rrt(032) * density(03) 
  pd(09,03) = pd(09,03) + rrt(032) * density(01) 
  pd(37,01) = pd(37,01) + rrt(032) * density(03) 
  pd(37,03) = pd(37,03) + rrt(032) * density(01) 
  pd(01,01) = pd(01,01) + rrt(033) * density(04) 
  pd(01,04) = pd(01,04) + rrt(033) * density(01) 
  pd(04,01) = pd(04,01) - rrt(033) * density(04) 
  pd(04,04) = pd(04,04) - rrt(033) * density(01) 
  pd(09,01) = pd(09,01) + rrt(033) * density(04) 
  pd(09,04) = pd(09,04) + rrt(033) * density(01) 
  pd(37,01) = pd(37,01) + rrt(033) * density(04) 
  pd(37,04) = pd(37,04) + rrt(033) * density(01) 
  pd(01,01) = pd(01,01) + rrt(034) * density(06) 
  pd(01,06) = pd(01,06) + rrt(034) * density(01) 
  pd(06,01) = pd(06,01) - rrt(034) * density(06) 
  pd(06,06) = pd(06,06) - rrt(034) * density(01) 
  pd(09,01) = pd(09,01) + rrt(034) * density(06) 
  pd(09,06) = pd(09,06) + rrt(034) * density(01) 
  pd(37,01) = pd(37,01) + rrt(034) * density(06) 
  pd(37,06) = pd(37,06) + rrt(034) * density(01) 
  pd(01,01) = pd(01,01) + rrt(035) * density(05) 
  pd(01,05) = pd(01,05) + rrt(035) * density(01) 
  pd(05,01) = pd(05,01) - rrt(035) * density(05) 
  pd(05,05) = pd(05,05) - rrt(035) * density(01) 
  pd(09,01) = pd(09,01) + rrt(035) * density(05) 
  pd(09,05) = pd(09,05) + rrt(035) * density(01) 
  pd(37,01) = pd(37,01) + rrt(035) * density(05) 
  pd(37,05) = pd(37,05) + rrt(035) * density(01) 
  pd(01,01) = pd(01,01) - rrt(036) * density(03) 
  pd(01,03) = pd(01,03) - rrt(036) * density(01) 
  pd(03,01) = pd(03,01) - rrt(036) * density(03) 
  pd(03,03) = pd(03,03) - rrt(036) * density(01) 
  pd(09,01) = pd(09,01) + rrt(036) * density(03) 
  pd(09,03) = pd(09,03) + rrt(036) * density(01) 
  pd(39,01) = pd(39,01) + rrt(036) * density(03) 
  pd(39,03) = pd(39,03) + rrt(036) * density(01) 
  pd(01,01) = pd(01,01) - rrt(037) * density(04) 
  pd(01,04) = pd(01,04) - rrt(037) * density(01) 
  pd(04,01) = pd(04,01) - rrt(037) * density(04) 
  pd(04,04) = pd(04,04) - rrt(037) * density(01) 
  pd(09,01) = pd(09,01) + rrt(037) * density(04) 
  pd(09,04) = pd(09,04) + rrt(037) * density(01) 
  pd(39,01) = pd(39,01) + rrt(037) * density(04) 
  pd(39,04) = pd(39,04) + rrt(037) * density(01) 
  pd(01,01) = pd(01,01) - rrt(038) * density(06) 
  pd(01,06) = pd(01,06) - rrt(038) * density(01) 
  pd(06,01) = pd(06,01) - rrt(038) * density(06) 
  pd(06,06) = pd(06,06) - rrt(038) * density(01) 
  pd(09,01) = pd(09,01) + rrt(038) * density(06) 
  pd(09,06) = pd(09,06) + rrt(038) * density(01) 
  pd(39,01) = pd(39,01) + rrt(038) * density(06) 
  pd(39,06) = pd(39,06) + rrt(038) * density(01) 
  pd(01,01) = pd(01,01) - rrt(039) * density(05) 
  pd(01,05) = pd(01,05) - rrt(039) * density(01) 
  pd(05,01) = pd(05,01) - rrt(039) * density(05) 
  pd(05,05) = pd(05,05) - rrt(039) * density(01) 
  pd(09,01) = pd(09,01) + rrt(039) * density(05) 
  pd(09,05) = pd(09,05) + rrt(039) * density(01) 
  pd(39,01) = pd(39,01) + rrt(039) * density(05) 
  pd(39,05) = pd(39,05) + rrt(039) * density(01) 
  pd(03,01) = pd(03,01) - rrt(040) * density(03) 
  pd(03,03) = pd(03,03) - rrt(040) * density(01) 
  pd(07,01) = pd(07,01) + rrt(040) * density(03) 
  pd(07,03) = pd(07,03) + rrt(040) * density(01) 
  pd(04,01) = pd(04,01) - rrt(041) * density(04) 
  pd(04,04) = pd(04,04) - rrt(041) * density(01) 
  pd(07,01) = pd(07,01) + rrt(041) * density(04) 
  pd(07,04) = pd(07,04) + rrt(041) * density(01) 
  pd(06,01) = pd(06,01) - rrt(042) * density(06) 
  pd(06,06) = pd(06,06) - rrt(042) * density(01) 
  pd(07,01) = pd(07,01) + rrt(042) * density(06) 
  pd(07,06) = pd(07,06) + rrt(042) * density(01) 
  pd(05,01) = pd(05,01) - rrt(043) * density(05) 
  pd(05,05) = pd(05,05) - rrt(043) * density(01) 
  pd(07,01) = pd(07,01) + rrt(043) * density(05) 
  pd(07,05) = pd(07,05) + rrt(043) * density(01) 
  pd(03,01) = pd(03,01) - rrt(044) * density(03) 
  pd(03,03) = pd(03,03) - rrt(044) * density(01) 
  pd(08,01) = pd(08,01) + rrt(044) * density(03) 
  pd(08,03) = pd(08,03) + rrt(044) * density(01) 
  pd(04,01) = pd(04,01) - rrt(045) * density(04) 
  pd(04,04) = pd(04,04) - rrt(045) * density(01) 
  pd(08,01) = pd(08,01) + rrt(045) * density(04) 
  pd(08,04) = pd(08,04) + rrt(045) * density(01) 
  pd(06,01) = pd(06,01) - rrt(046) * density(06) 
  pd(06,06) = pd(06,06) - rrt(046) * density(01) 
  pd(08,01) = pd(08,01) + rrt(046) * density(06) 
  pd(08,06) = pd(08,06) + rrt(046) * density(01) 
  pd(05,01) = pd(05,01) - rrt(047) * density(05) 
  pd(05,05) = pd(05,05) - rrt(047) * density(01) 
  pd(08,01) = pd(08,01) + rrt(047) * density(05) 
  pd(08,05) = pd(08,05) + rrt(047) * density(01) 
  pd(01,01) = pd(01,01) + rrt(048) * density(03) 
  pd(01,03) = pd(01,03) + rrt(048) * density(01) 
  pd(03,01) = pd(03,01) - rrt(048) * density(03) 
  pd(03,03) = pd(03,03) - rrt(048) * density(01) 
  pd(22,01) = pd(22,01) + rrt(048) * density(03) 
  pd(22,03) = pd(22,03) + rrt(048) * density(01) 
  pd(36,01) = pd(36,01) + rrt(048) * density(03) 
  pd(36,03) = pd(36,03) + rrt(048) * density(01) 
  pd(01,01) = pd(01,01) + rrt(049) * density(04) 
  pd(01,04) = pd(01,04) + rrt(049) * density(01) 
  pd(04,01) = pd(04,01) - rrt(049) * density(04) 
  pd(04,04) = pd(04,04) - rrt(049) * density(01) 
  pd(22,01) = pd(22,01) + rrt(049) * density(04) 
  pd(22,04) = pd(22,04) + rrt(049) * density(01) 
  pd(36,01) = pd(36,01) + rrt(049) * density(04) 
  pd(36,04) = pd(36,04) + rrt(049) * density(01) 
  pd(01,01) = pd(01,01) + rrt(050) * density(06) 
  pd(01,06) = pd(01,06) + rrt(050) * density(01) 
  pd(06,01) = pd(06,01) - rrt(050) * density(06) 
  pd(06,06) = pd(06,06) - rrt(050) * density(01) 
  pd(22,01) = pd(22,01) + rrt(050) * density(06) 
  pd(22,06) = pd(22,06) + rrt(050) * density(01) 
  pd(36,01) = pd(36,01) + rrt(050) * density(06) 
  pd(36,06) = pd(36,06) + rrt(050) * density(01) 
  pd(01,01) = pd(01,01) + rrt(051) * density(05) 
  pd(01,05) = pd(01,05) + rrt(051) * density(01) 
  pd(05,01) = pd(05,01) - rrt(051) * density(05) 
  pd(05,05) = pd(05,05) - rrt(051) * density(01) 
  pd(22,01) = pd(22,01) + rrt(051) * density(05) 
  pd(22,05) = pd(22,05) + rrt(051) * density(01) 
  pd(36,01) = pd(36,01) + rrt(051) * density(05) 
  pd(36,05) = pd(36,05) + rrt(051) * density(01) 
  pd(02,01) = pd(02,01) + rrt(052) * density(03) 
  pd(02,03) = pd(02,03) + rrt(052) * density(01) 
  pd(03,01) = pd(03,01) - rrt(052) * density(03) 
  pd(03,03) = pd(03,03) - rrt(052) * density(01) 
  pd(02,01) = pd(02,01) + rrt(053) * density(04) 
  pd(02,04) = pd(02,04) + rrt(053) * density(01) 
  pd(04,01) = pd(04,01) - rrt(053) * density(04) 
  pd(04,04) = pd(04,04) - rrt(053) * density(01) 
  pd(02,01) = pd(02,01) + rrt(054) * density(06) 
  pd(02,06) = pd(02,06) + rrt(054) * density(01) 
  pd(06,01) = pd(06,01) - rrt(054) * density(06) 
  pd(06,06) = pd(06,06) - rrt(054) * density(01) 
  pd(02,01) = pd(02,01) + rrt(055) * density(05) 
  pd(02,05) = pd(02,05) + rrt(055) * density(01) 
  pd(05,01) = pd(05,01) - rrt(055) * density(05) 
  pd(05,05) = pd(05,05) - rrt(055) * density(01) 
  pd(07,01) = pd(07,01) - rrt(056) * density(07) 
  pd(07,07) = pd(07,07) - rrt(056) * density(01) 
  pd(09,01) = pd(09,01) + rrt(056) * density(07) 
  pd(09,07) = pd(09,07) + rrt(056) * density(01) 
  pd(25,01) = pd(25,01) + rrt(056) * density(07) 
  pd(25,07) = pd(25,07) + rrt(056) * density(01) 
  pd(08,01) = pd(08,01) - rrt(057) * density(08) 
  pd(08,08) = pd(08,08) - rrt(057) * density(01) 
  pd(09,01) = pd(09,01) + rrt(057) * density(08) 
  pd(09,08) = pd(09,08) + rrt(057) * density(01) 
  pd(25,01) = pd(25,01) + rrt(057) * density(08) 
  pd(25,08) = pd(25,08) + rrt(057) * density(01) 
  pd(01,01) = pd(01,01) + rrt(058) * density(07) 
  pd(01,07) = pd(01,07) + rrt(058) * density(01) 
  pd(07,01) = pd(07,01) - rrt(058) * density(07) 
  pd(07,07) = pd(07,07) - rrt(058) * density(01) 
  pd(26,01) = pd(26,01) + rrt(058) * density(07) 
  pd(26,07) = pd(26,07) + rrt(058) * density(01) 
  pd(01,01) = pd(01,01) + rrt(059) * density(08) 
  pd(01,08) = pd(01,08) + rrt(059) * density(01) 
  pd(08,01) = pd(08,01) - rrt(059) * density(08) 
  pd(08,08) = pd(08,08) - rrt(059) * density(01) 
  pd(26,01) = pd(26,01) + rrt(059) * density(08) 
  pd(26,08) = pd(26,08) + rrt(059) * density(01) 
  pd(01,01) = pd(01,01) + rrt(060) * density(07) 
  pd(01,07) = pd(01,07) + rrt(060) * density(01) 
  pd(07,01) = pd(07,01) - rrt(060) * density(07) 
  pd(07,07) = pd(07,07) - rrt(060) * density(01) 
  pd(25,01) = pd(25,01) + rrt(060) * density(07) 
  pd(25,07) = pd(25,07) + rrt(060) * density(01) 
  pd(28,01) = pd(28,01) + rrt(060) * density(07) 
  pd(28,07) = pd(28,07) + rrt(060) * density(01) 
  pd(01,01) = pd(01,01) + rrt(061) * density(08) 
  pd(01,08) = pd(01,08) + rrt(061) * density(01) 
  pd(08,01) = pd(08,01) - rrt(061) * density(08) 
  pd(08,08) = pd(08,08) - rrt(061) * density(01) 
  pd(25,01) = pd(25,01) + rrt(061) * density(08) 
  pd(25,08) = pd(25,08) + rrt(061) * density(01) 
  pd(28,01) = pd(28,01) + rrt(061) * density(08) 
  pd(28,08) = pd(28,08) + rrt(061) * density(01) 
  pd(01,01) = pd(01,01) + rrt(062) * density(07) 
  pd(01,07) = pd(01,07) + rrt(062) * density(01) 
  pd(07,01) = pd(07,01) - rrt(062) * density(07) 
  pd(07,07) = pd(07,07) - rrt(062) * density(01) 
  pd(15,01) = pd(15,01) + rrt(062) * density(07) 
  pd(15,07) = pd(15,07) + rrt(062) * density(01) 
  pd(32,01) = pd(32,01) + rrt(062) * density(07) 
  pd(32,07) = pd(32,07) + rrt(062) * density(01) 
  pd(01,01) = pd(01,01) + rrt(063) * density(08) 
  pd(01,08) = pd(01,08) + rrt(063) * density(01) 
  pd(08,01) = pd(08,01) - rrt(063) * density(08) 
  pd(08,08) = pd(08,08) - rrt(063) * density(01) 
  pd(15,01) = pd(15,01) + rrt(063) * density(08) 
  pd(15,08) = pd(15,08) + rrt(063) * density(01) 
  pd(32,01) = pd(32,01) + rrt(063) * density(08) 
  pd(32,08) = pd(32,08) + rrt(063) * density(01) 
  pd(01,01) = pd(01,01) + rrt(064) * density(07) 
  pd(01,07) = pd(01,07) + rrt(064) * density(01) 
  pd(07,01) = pd(07,01) - rrt(064) * density(07) 
  pd(07,07) = pd(07,07) - rrt(064) * density(01) 
  pd(09,01) = pd(09,01) + rrt(064) * density(07) 
  pd(09,07) = pd(09,07) + rrt(064) * density(01) 
  pd(37,01) = pd(37,01) + rrt(064) * density(07) 
  pd(37,07) = pd(37,07) + rrt(064) * density(01) 
  pd(01,01) = pd(01,01) + rrt(065) * density(08) 
  pd(01,08) = pd(01,08) + rrt(065) * density(01) 
  pd(08,01) = pd(08,01) - rrt(065) * density(08) 
  pd(08,08) = pd(08,08) - rrt(065) * density(01) 
  pd(09,01) = pd(09,01) + rrt(065) * density(08) 
  pd(09,08) = pd(09,08) + rrt(065) * density(01) 
  pd(37,01) = pd(37,01) + rrt(065) * density(08) 
  pd(37,08) = pd(37,08) + rrt(065) * density(01) 
  pd(01,01) = pd(01,01) + rrt(066) * density(07) 
  pd(01,07) = pd(01,07) + rrt(066) * density(01) 
  pd(07,01) = pd(07,01) - rrt(066) * density(07) 
  pd(07,07) = pd(07,07) - rrt(066) * density(01) 
  pd(22,01) = pd(22,01) + rrt(066) * density(07) 
  pd(22,07) = pd(22,07) + rrt(066) * density(01) 
  pd(36,01) = pd(36,01) + rrt(066) * density(07) 
  pd(36,07) = pd(36,07) + rrt(066) * density(01) 
  pd(01,01) = pd(01,01) + rrt(067) * density(08) 
  pd(01,08) = pd(01,08) + rrt(067) * density(01) 
  pd(08,01) = pd(08,01) - rrt(067) * density(08) 
  pd(08,08) = pd(08,08) - rrt(067) * density(01) 
  pd(22,01) = pd(22,01) + rrt(067) * density(08) 
  pd(22,08) = pd(22,08) + rrt(067) * density(01) 
  pd(36,01) = pd(36,01) + rrt(067) * density(08) 
  pd(36,08) = pd(36,08) + rrt(067) * density(01) 
  pd(01,01) = pd(01,01) + rrt(068) * density(09) 
  pd(01,09) = pd(01,09) + rrt(068) * density(01) 
  pd(09,01) = pd(09,01) - rrt(068) * density(09) 
  pd(09,09) = pd(09,09) - rrt(068) * density(01) 
  pd(28,01) = pd(28,01) + rrt(068) * density(09) 
  pd(28,09) = pd(28,09) + rrt(068) * density(01) 
  pd(01,01) = pd(01,01) + rrt(069) * density(09) 
  pd(01,09) = pd(01,09) + rrt(069) * density(01) 
  pd(09,01) = pd(09,01) - rrt(069) * density(09) 
  pd(09,09) = pd(09,09) - rrt(069) * density(01) 
  pd(25,01) = pd(25,01) + rrt(069) * density(09) 
  pd(25,09) = pd(25,09) + rrt(069) * density(01) 
  pd(32,01) = pd(32,01) + rrt(069) * density(09) 
  pd(32,09) = pd(32,09) + rrt(069) * density(01) 
  pd(01,01) = pd(01,01) + rrt(070) * density(09) 
  pd(01,09) = pd(01,09) + rrt(070) * density(01) 
  pd(09,01) = pd(09,01) - rrt(070) * density(09) 
  pd(09,09) = pd(09,09) - rrt(070) * density(01) 
  pd(22,01) = pd(22,01) + rrt(070) * density(09) 
  pd(22,09) = pd(22,09) + rrt(070) * density(01) 
  pd(37,01) = pd(37,01) + rrt(070) * density(09) 
  pd(37,09) = pd(37,09) + rrt(070) * density(01) 
  pd(01,01) = pd(01,01) - rrt(071) * density(09) 
  pd(01,09) = pd(01,09) - rrt(071) * density(01) 
  pd(09,01) = pd(09,01) - rrt(071) * density(09) 
  pd(09,09) = pd(09,09) - rrt(071) * density(01) 
  pd(22,01) = pd(22,01) + rrt(071) * density(09) 
  pd(22,09) = pd(22,09) + rrt(071) * density(01) 
  pd(39,01) = pd(39,01) + rrt(071) * density(09) 
  pd(39,09) = pd(39,09) + rrt(071) * density(01) 
  pd(09,01) = pd(09,01) - rrt(072) * density(09) 
  pd(09,09) = pd(09,09) - rrt(072) * density(01) 
  pd(22,01) = pd(22,01) + rrt(072) * density(09) 
  pd(22,09) = pd(22,09) + rrt(072) * density(01) 
  pd(25,01) = pd(25,01) + rrt(072) * density(09) 
  pd(25,09) = pd(25,09) + rrt(072) * density(01) 
  pd(09,01) = pd(09,01) - rrt(073) * density(09) 
  pd(09,09) = pd(09,09) - rrt(073) * density(01) 
  pd(10,01) = pd(10,01) + rrt(073) * density(09) 
  pd(10,09) = pd(10,09) + rrt(073) * density(01) 
  pd(09,01) = pd(09,01) - rrt(074) * density(09) 
  pd(09,09) = pd(09,09) - rrt(074) * density(01) 
  pd(11,01) = pd(11,01) + rrt(074) * density(09) 
  pd(11,09) = pd(11,09) + rrt(074) * density(01) 
  pd(09,01) = pd(09,01) - rrt(075) * density(09) 
  pd(09,09) = pd(09,09) - rrt(075) * density(01) 
  pd(12,01) = pd(12,01) + rrt(075) * density(09) 
  pd(12,09) = pd(12,09) + rrt(075) * density(01) 
  pd(09,01) = pd(09,01) - rrt(076) * density(09) 
  pd(09,09) = pd(09,09) - rrt(076) * density(01) 
  pd(13,01) = pd(13,01) + rrt(076) * density(09) 
  pd(13,09) = pd(13,09) + rrt(076) * density(01) 
  pd(09,01) = pd(09,01) - rrt(077) * density(09) 
  pd(09,09) = pd(09,09) - rrt(077) * density(01) 
  pd(14,01) = pd(14,01) + rrt(077) * density(09) 
  pd(14,09) = pd(14,09) + rrt(077) * density(01) 
  pd(09,01) = pd(09,01) + rrt(078) * density(11) 
  pd(09,11) = pd(09,11) + rrt(078) * density(01) 
  pd(11,01) = pd(11,01) - rrt(078) * density(11) 
  pd(11,11) = pd(11,11) - rrt(078) * density(01) 
  pd(09,01) = pd(09,01) + rrt(079) * density(12) 
  pd(09,12) = pd(09,12) + rrt(079) * density(01) 
  pd(12,01) = pd(12,01) - rrt(079) * density(12) 
  pd(12,12) = pd(12,12) - rrt(079) * density(01) 
  pd(09,01) = pd(09,01) + rrt(080) * density(13) 
  pd(09,13) = pd(09,13) + rrt(080) * density(01) 
  pd(13,01) = pd(13,01) - rrt(080) * density(13) 
  pd(13,13) = pd(13,13) - rrt(080) * density(01) 
  pd(09,01) = pd(09,01) + rrt(081) * density(14) 
  pd(09,14) = pd(09,14) + rrt(081) * density(01) 
  pd(14,01) = pd(14,01) - rrt(081) * density(14) 
  pd(14,14) = pd(14,14) - rrt(081) * density(01) 
  pd(10,01) = pd(10,01) - rrt(082) * density(10) 
  pd(10,10) = pd(10,10) - rrt(082) * density(01) 
  pd(22,01) = pd(22,01) + rrt(082) * density(10) 
  pd(22,10) = pd(22,10) + rrt(082) * density(01) 
  pd(25,01) = pd(25,01) + rrt(082) * density(10) 
  pd(25,10) = pd(25,10) + rrt(082) * density(01) 
  pd(01,01) = pd(01,01) + rrt(083) * density(10) 
  pd(01,10) = pd(01,10) + rrt(083) * density(01) 
  pd(10,01) = pd(10,01) - rrt(083) * density(10) 
  pd(10,10) = pd(10,10) - rrt(083) * density(01) 
  pd(28,01) = pd(28,01) + rrt(083) * density(10) 
  pd(28,10) = pd(28,10) + rrt(083) * density(01) 
  pd(01,01) = pd(01,01) + rrt(084) * density(10) 
  pd(01,10) = pd(01,10) + rrt(084) * density(01) 
  pd(10,01) = pd(10,01) - rrt(084) * density(10) 
  pd(10,10) = pd(10,10) - rrt(084) * density(01) 
  pd(25,01) = pd(25,01) + rrt(084) * density(10) 
  pd(25,10) = pd(25,10) + rrt(084) * density(01) 
  pd(32,01) = pd(32,01) + rrt(084) * density(10) 
  pd(32,10) = pd(32,10) + rrt(084) * density(01) 
  pd(01,01) = pd(01,01) + rrt(085) * density(10) 
  pd(01,10) = pd(01,10) + rrt(085) * density(01) 
  pd(10,01) = pd(10,01) - rrt(085) * density(10) 
  pd(10,10) = pd(10,10) - rrt(085) * density(01) 
  pd(22,01) = pd(22,01) + rrt(085) * density(10) 
  pd(22,10) = pd(22,10) + rrt(085) * density(01) 
  pd(37,01) = pd(37,01) + rrt(085) * density(10) 
  pd(37,10) = pd(37,10) + rrt(085) * density(01) 
  pd(01,01) = pd(01,01) - rrt(086) * density(10) 
  pd(01,10) = pd(01,10) - rrt(086) * density(01) 
  pd(10,01) = pd(10,01) - rrt(086) * density(10) 
  pd(10,10) = pd(10,10) - rrt(086) * density(01) 
  pd(22,01) = pd(22,01) + rrt(086) * density(10) 
  pd(22,10) = pd(22,10) + rrt(086) * density(01) 
  pd(39,01) = pd(39,01) + rrt(086) * density(10) 
  pd(39,10) = pd(39,10) + rrt(086) * density(01) 
  pd(10,01) = pd(10,01) - rrt(087) * density(10) 
  pd(10,10) = pd(10,10) - rrt(087) * density(01) 
  pd(11,01) = pd(11,01) + rrt(087) * density(10) 
  pd(11,10) = pd(11,10) + rrt(087) * density(01) 
  pd(10,01) = pd(10,01) - rrt(088) * density(10) 
  pd(10,10) = pd(10,10) - rrt(088) * density(01) 
  pd(12,01) = pd(12,01) + rrt(088) * density(10) 
  pd(12,10) = pd(12,10) + rrt(088) * density(01) 
  pd(10,01) = pd(10,01) - rrt(089) * density(10) 
  pd(10,10) = pd(10,10) - rrt(089) * density(01) 
  pd(13,01) = pd(13,01) + rrt(089) * density(10) 
  pd(13,10) = pd(13,10) + rrt(089) * density(01) 
  pd(10,01) = pd(10,01) - rrt(090) * density(10) 
  pd(10,10) = pd(10,10) - rrt(090) * density(01) 
  pd(14,01) = pd(14,01) + rrt(090) * density(10) 
  pd(14,10) = pd(14,10) + rrt(090) * density(01) 
  pd(09,01) = pd(09,01) + rrt(091) * density(10) 
  pd(09,10) = pd(09,10) + rrt(091) * density(01) 
  pd(10,01) = pd(10,01) - rrt(091) * density(10) 
  pd(10,10) = pd(10,10) - rrt(091) * density(01) 
  pd(11,01) = pd(11,01) - rrt(092) * density(11) 
  pd(11,11) = pd(11,11) - rrt(092) * density(01) 
  pd(22,01) = pd(22,01) + rrt(092) * density(11) 
  pd(22,11) = pd(22,11) + rrt(092) * density(01) 
  pd(25,01) = pd(25,01) + rrt(092) * density(11) 
  pd(25,11) = pd(25,11) + rrt(092) * density(01) 
  pd(12,01) = pd(12,01) - rrt(093) * density(12) 
  pd(12,12) = pd(12,12) - rrt(093) * density(01) 
  pd(22,01) = pd(22,01) + rrt(093) * density(12) 
  pd(22,12) = pd(22,12) + rrt(093) * density(01) 
  pd(25,01) = pd(25,01) + rrt(093) * density(12) 
  pd(25,12) = pd(25,12) + rrt(093) * density(01) 
  pd(13,01) = pd(13,01) - rrt(094) * density(13) 
  pd(13,13) = pd(13,13) - rrt(094) * density(01) 
  pd(22,01) = pd(22,01) + rrt(094) * density(13) 
  pd(22,13) = pd(22,13) + rrt(094) * density(01) 
  pd(25,01) = pd(25,01) + rrt(094) * density(13) 
  pd(25,13) = pd(25,13) + rrt(094) * density(01) 
  pd(14,01) = pd(14,01) - rrt(095) * density(14) 
  pd(14,14) = pd(14,14) - rrt(095) * density(01) 
  pd(22,01) = pd(22,01) + rrt(095) * density(14) 
  pd(22,14) = pd(22,14) + rrt(095) * density(01) 
  pd(25,01) = pd(25,01) + rrt(095) * density(14) 
  pd(25,14) = pd(25,14) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(11) 
  pd(01,11) = pd(01,11) + rrt(096) * density(01) 
  pd(11,01) = pd(11,01) - rrt(096) * density(11) 
  pd(11,11) = pd(11,11) - rrt(096) * density(01) 
  pd(28,01) = pd(28,01) + rrt(096) * density(11) 
  pd(28,11) = pd(28,11) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(12) 
  pd(01,12) = pd(01,12) + rrt(097) * density(01) 
  pd(12,01) = pd(12,01) - rrt(097) * density(12) 
  pd(12,12) = pd(12,12) - rrt(097) * density(01) 
  pd(28,01) = pd(28,01) + rrt(097) * density(12) 
  pd(28,12) = pd(28,12) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(13) 
  pd(01,13) = pd(01,13) + rrt(098) * density(01) 
  pd(13,01) = pd(13,01) - rrt(098) * density(13) 
  pd(13,13) = pd(13,13) - rrt(098) * density(01) 
  pd(28,01) = pd(28,01) + rrt(098) * density(13) 
  pd(28,13) = pd(28,13) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(14) 
  pd(01,14) = pd(01,14) + rrt(099) * density(01) 
  pd(14,01) = pd(14,01) - rrt(099) * density(14) 
  pd(14,14) = pd(14,14) - rrt(099) * density(01) 
  pd(28,01) = pd(28,01) + rrt(099) * density(14) 
  pd(28,14) = pd(28,14) + rrt(099) * density(01) 
  pd(01,01) = pd(01,01) + rrt(100) * density(11) 
  pd(01,11) = pd(01,11) + rrt(100) * density(01) 
  pd(11,01) = pd(11,01) - rrt(100) * density(11) 
  pd(11,11) = pd(11,11) - rrt(100) * density(01) 
  pd(25,01) = pd(25,01) + rrt(100) * density(11) 
  pd(25,11) = pd(25,11) + rrt(100) * density(01) 
  pd(32,01) = pd(32,01) + rrt(100) * density(11) 
  pd(32,11) = pd(32,11) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(12) 
  pd(01,12) = pd(01,12) + rrt(101) * density(01) 
  pd(12,01) = pd(12,01) - rrt(101) * density(12) 
  pd(12,12) = pd(12,12) - rrt(101) * density(01) 
  pd(25,01) = pd(25,01) + rrt(101) * density(12) 
  pd(25,12) = pd(25,12) + rrt(101) * density(01) 
  pd(32,01) = pd(32,01) + rrt(101) * density(12) 
  pd(32,12) = pd(32,12) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(13) 
  pd(01,13) = pd(01,13) + rrt(102) * density(01) 
  pd(13,01) = pd(13,01) - rrt(102) * density(13) 
  pd(13,13) = pd(13,13) - rrt(102) * density(01) 
  pd(25,01) = pd(25,01) + rrt(102) * density(13) 
  pd(25,13) = pd(25,13) + rrt(102) * density(01) 
  pd(32,01) = pd(32,01) + rrt(102) * density(13) 
  pd(32,13) = pd(32,13) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(14) 
  pd(01,14) = pd(01,14) + rrt(103) * density(01) 
  pd(14,01) = pd(14,01) - rrt(103) * density(14) 
  pd(14,14) = pd(14,14) - rrt(103) * density(01) 
  pd(25,01) = pd(25,01) + rrt(103) * density(14) 
  pd(25,14) = pd(25,14) + rrt(103) * density(01) 
  pd(32,01) = pd(32,01) + rrt(103) * density(14) 
  pd(32,14) = pd(32,14) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(11) 
  pd(01,11) = pd(01,11) + rrt(104) * density(01) 
  pd(11,01) = pd(11,01) - rrt(104) * density(11) 
  pd(11,11) = pd(11,11) - rrt(104) * density(01) 
  pd(22,01) = pd(22,01) + rrt(104) * density(11) 
  pd(22,11) = pd(22,11) + rrt(104) * density(01) 
  pd(37,01) = pd(37,01) + rrt(104) * density(11) 
  pd(37,11) = pd(37,11) + rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(12) 
  pd(01,12) = pd(01,12) + rrt(105) * density(01) 
  pd(12,01) = pd(12,01) - rrt(105) * density(12) 
  pd(12,12) = pd(12,12) - rrt(105) * density(01) 
  pd(22,01) = pd(22,01) + rrt(105) * density(12) 
  pd(22,12) = pd(22,12) + rrt(105) * density(01) 
  pd(37,01) = pd(37,01) + rrt(105) * density(12) 
  pd(37,12) = pd(37,12) + rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(13) 
  pd(01,13) = pd(01,13) + rrt(106) * density(01) 
  pd(13,01) = pd(13,01) - rrt(106) * density(13) 
  pd(13,13) = pd(13,13) - rrt(106) * density(01) 
  pd(22,01) = pd(22,01) + rrt(106) * density(13) 
  pd(22,13) = pd(22,13) + rrt(106) * density(01) 
  pd(37,01) = pd(37,01) + rrt(106) * density(13) 
  pd(37,13) = pd(37,13) + rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(14) 
  pd(01,14) = pd(01,14) + rrt(107) * density(01) 
  pd(14,01) = pd(14,01) - rrt(107) * density(14) 
  pd(14,14) = pd(14,14) - rrt(107) * density(01) 
  pd(22,01) = pd(22,01) + rrt(107) * density(14) 
  pd(22,14) = pd(22,14) + rrt(107) * density(01) 
  pd(37,01) = pd(37,01) + rrt(107) * density(14) 
  pd(37,14) = pd(37,14) + rrt(107) * density(01) 
  pd(01,01) = pd(01,01) + rrt(108) * density(22) 
  pd(01,22) = pd(01,22) + rrt(108) * density(01) 
  pd(22,01) = pd(22,01) - rrt(108) * density(22) 
  pd(22,22) = pd(22,22) - rrt(108) * density(01) 
  pd(32,01) = pd(32,01) + rrt(108) * density(22) 
  pd(32,22) = pd(32,22) + rrt(108) * density(01) 
  pd(22,01) = pd(22,01) + rrt(109) * density(23) * 2.0d0
  pd(22,23) = pd(22,23) + rrt(109) * density(01) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(109) * density(23) 
  pd(23,23) = pd(23,23) - rrt(109) * density(01) 
  pd(01,01) = pd(01,01) + rrt(110) * density(23) 
  pd(01,23) = pd(01,23) + rrt(110) * density(01) 
  pd(23,01) = pd(23,01) - rrt(110) * density(23) 
  pd(23,23) = pd(23,23) - rrt(110) * density(01) 
  pd(33,01) = pd(33,01) + rrt(110) * density(23) 
  pd(33,23) = pd(33,23) + rrt(110) * density(01) 
  pd(15,01) = pd(15,01) - rrt(111) * density(15) 
  pd(15,15) = pd(15,15) - rrt(111) * density(01) 
  pd(25,01) = pd(25,01) + rrt(111) * density(15) * 2.0d0
  pd(25,15) = pd(25,15) + rrt(111) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(112) * density(15) 
  pd(01,15) = pd(01,15) + rrt(112) * density(01) 
  pd(15,01) = pd(15,01) - rrt(112) * density(15) 
  pd(15,15) = pd(15,15) - rrt(112) * density(01) 
  pd(36,01) = pd(36,01) + rrt(112) * density(15) 
  pd(36,15) = pd(36,15) + rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(15) 
  pd(01,15) = pd(01,15) + rrt(113) * density(01) 
  pd(15,01) = pd(15,01) - rrt(113) * density(15) 
  pd(15,15) = pd(15,15) - rrt(113) * density(01) 
  pd(25,01) = pd(25,01) + rrt(113) * density(15) 
  pd(25,15) = pd(25,15) + rrt(113) * density(01) 
  pd(37,01) = pd(37,01) + rrt(113) * density(15) 
  pd(37,15) = pd(37,15) + rrt(113) * density(01) 
  pd(01,01) = pd(01,01) - rrt(114) * density(15) 
  pd(01,15) = pd(01,15) - rrt(114) * density(01) 
  pd(15,01) = pd(15,01) - rrt(114) * density(15) 
  pd(15,15) = pd(15,15) - rrt(114) * density(01) 
  pd(25,01) = pd(25,01) + rrt(114) * density(15) 
  pd(25,15) = pd(25,15) + rrt(114) * density(01) 
  pd(39,01) = pd(39,01) + rrt(114) * density(15) 
  pd(39,15) = pd(39,15) + rrt(114) * density(01) 
  pd(15,01) = pd(15,01) - rrt(115) * density(15) 
  pd(15,15) = pd(15,15) - rrt(115) * density(01) 
  pd(16,01) = pd(16,01) + rrt(115) * density(15) 
  pd(16,15) = pd(16,15) + rrt(115) * density(01) 
  pd(15,01) = pd(15,01) - rrt(116) * density(15) 
  pd(15,15) = pd(15,15) - rrt(116) * density(01) 
  pd(17,01) = pd(17,01) + rrt(116) * density(15) 
  pd(17,15) = pd(17,15) + rrt(116) * density(01) 
  pd(15,01) = pd(15,01) - rrt(117) * density(15) 
  pd(15,15) = pd(15,15) - rrt(117) * density(01) 
  pd(18,01) = pd(18,01) + rrt(117) * density(15) 
  pd(18,15) = pd(18,15) + rrt(117) * density(01) 
  pd(15,01) = pd(15,01) - rrt(118) * density(15) 
  pd(15,15) = pd(15,15) - rrt(118) * density(01) 
  pd(19,01) = pd(19,01) + rrt(118) * density(15) 
  pd(19,15) = pd(19,15) + rrt(118) * density(01) 
  pd(15,01) = pd(15,01) - rrt(119) * density(15) 
  pd(15,15) = pd(15,15) - rrt(119) * density(01) 
  pd(20,01) = pd(20,01) + rrt(119) * density(15) 
  pd(20,15) = pd(20,15) + rrt(119) * density(01) 
  pd(15,01) = pd(15,01) + rrt(120) * density(19) 
  pd(15,19) = pd(15,19) + rrt(120) * density(01) 
  pd(19,01) = pd(19,01) - rrt(120) * density(19) 
  pd(19,19) = pd(19,19) - rrt(120) * density(01) 
  pd(15,01) = pd(15,01) + rrt(121) * density(20) 
  pd(15,20) = pd(15,20) + rrt(121) * density(01) 
  pd(20,01) = pd(20,01) - rrt(121) * density(20) 
  pd(20,20) = pd(20,20) - rrt(121) * density(01) 
  pd(15,19) = pd(15,19) + rrt(122) 
  pd(19,19) = pd(19,19) - rrt(122) 
  pd(15,20) = pd(15,20) + rrt(123) 
  pd(20,20) = pd(20,20) - rrt(123) 
  pd(16,01) = pd(16,01) - rrt(124) * density(16) 
  pd(16,16) = pd(16,16) - rrt(124) * density(01) 
  pd(19,01) = pd(19,01) + rrt(124) * density(16) 
  pd(19,16) = pd(19,16) + rrt(124) * density(01) 
  pd(17,01) = pd(17,01) - rrt(125) * density(17) 
  pd(17,17) = pd(17,17) - rrt(125) * density(01) 
  pd(19,01) = pd(19,01) + rrt(125) * density(17) 
  pd(19,17) = pd(19,17) + rrt(125) * density(01) 
  pd(18,01) = pd(18,01) - rrt(126) * density(18) 
  pd(18,18) = pd(18,18) - rrt(126) * density(01) 
  pd(19,01) = pd(19,01) + rrt(126) * density(18) 
  pd(19,18) = pd(19,18) + rrt(126) * density(01) 
  pd(16,01) = pd(16,01) - rrt(127) * density(16) 
  pd(16,16) = pd(16,16) - rrt(127) * density(01) 
  pd(20,01) = pd(20,01) + rrt(127) * density(16) 
  pd(20,16) = pd(20,16) + rrt(127) * density(01) 
  pd(17,01) = pd(17,01) - rrt(128) * density(17) 
  pd(17,17) = pd(17,17) - rrt(128) * density(01) 
  pd(20,01) = pd(20,01) + rrt(128) * density(17) 
  pd(20,17) = pd(20,17) + rrt(128) * density(01) 
  pd(18,01) = pd(18,01) - rrt(129) * density(18) 
  pd(18,18) = pd(18,18) - rrt(129) * density(01) 
  pd(20,01) = pd(20,01) + rrt(129) * density(18) 
  pd(20,18) = pd(20,18) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(16) 
  pd(01,16) = pd(01,16) + rrt(130) * density(01) 
  pd(16,01) = pd(16,01) - rrt(130) * density(16) 
  pd(16,16) = pd(16,16) - rrt(130) * density(01) 
  pd(36,01) = pd(36,01) + rrt(130) * density(16) 
  pd(36,16) = pd(36,16) + rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(17) 
  pd(01,17) = pd(01,17) + rrt(131) * density(01) 
  pd(17,01) = pd(17,01) - rrt(131) * density(17) 
  pd(17,17) = pd(17,17) - rrt(131) * density(01) 
  pd(36,01) = pd(36,01) + rrt(131) * density(17) 
  pd(36,17) = pd(36,17) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(18) 
  pd(01,18) = pd(01,18) + rrt(132) * density(01) 
  pd(18,01) = pd(18,01) - rrt(132) * density(18) 
  pd(18,18) = pd(18,18) - rrt(132) * density(01) 
  pd(36,01) = pd(36,01) + rrt(132) * density(18) 
  pd(36,18) = pd(36,18) + rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(16) 
  pd(01,16) = pd(01,16) + rrt(133) * density(01) 
  pd(16,01) = pd(16,01) - rrt(133) * density(16) 
  pd(16,16) = pd(16,16) - rrt(133) * density(01) 
  pd(25,01) = pd(25,01) + rrt(133) * density(16) 
  pd(25,16) = pd(25,16) + rrt(133) * density(01) 
  pd(37,01) = pd(37,01) + rrt(133) * density(16) 
  pd(37,16) = pd(37,16) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(17) 
  pd(01,17) = pd(01,17) + rrt(134) * density(01) 
  pd(17,01) = pd(17,01) - rrt(134) * density(17) 
  pd(17,17) = pd(17,17) - rrt(134) * density(01) 
  pd(25,01) = pd(25,01) + rrt(134) * density(17) 
  pd(25,17) = pd(25,17) + rrt(134) * density(01) 
  pd(37,01) = pd(37,01) + rrt(134) * density(17) 
  pd(37,17) = pd(37,17) + rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(18) 
  pd(01,18) = pd(01,18) + rrt(135) * density(01) 
  pd(18,01) = pd(18,01) - rrt(135) * density(18) 
  pd(18,18) = pd(18,18) - rrt(135) * density(01) 
  pd(25,01) = pd(25,01) + rrt(135) * density(18) 
  pd(25,18) = pd(25,18) + rrt(135) * density(01) 
  pd(37,01) = pd(37,01) + rrt(135) * density(18) 
  pd(37,18) = pd(37,18) + rrt(135) * density(01) 
  pd(16,01) = pd(16,01) - rrt(136) * density(16) 
  pd(16,16) = pd(16,16) - rrt(136) * density(01) 
  pd(25,01) = pd(25,01) + rrt(136) * density(16) * 2.0d0
  pd(25,16) = pd(25,16) + rrt(136) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) - rrt(137) * density(17) 
  pd(17,17) = pd(17,17) - rrt(137) * density(01) 
  pd(25,01) = pd(25,01) + rrt(137) * density(17) * 2.0d0
  pd(25,17) = pd(25,17) + rrt(137) * density(01) * 2.0d0
  pd(18,01) = pd(18,01) - rrt(138) * density(18) 
  pd(18,18) = pd(18,18) - rrt(138) * density(01) 
  pd(25,01) = pd(25,01) + rrt(138) * density(18) * 2.0d0
  pd(25,18) = pd(25,18) + rrt(138) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(139) * density(16) 
  pd(01,16) = pd(01,16) - rrt(139) * density(01) 
  pd(16,01) = pd(16,01) - rrt(139) * density(16) 
  pd(16,16) = pd(16,16) - rrt(139) * density(01) 
  pd(25,01) = pd(25,01) + rrt(139) * density(16) 
  pd(25,16) = pd(25,16) + rrt(139) * density(01) 
  pd(39,01) = pd(39,01) + rrt(139) * density(16) 
  pd(39,16) = pd(39,16) + rrt(139) * density(01) 
  pd(01,01) = pd(01,01) - rrt(140) * density(17) 
  pd(01,17) = pd(01,17) - rrt(140) * density(01) 
  pd(17,01) = pd(17,01) - rrt(140) * density(17) 
  pd(17,17) = pd(17,17) - rrt(140) * density(01) 
  pd(25,01) = pd(25,01) + rrt(140) * density(17) 
  pd(25,17) = pd(25,17) + rrt(140) * density(01) 
  pd(39,01) = pd(39,01) + rrt(140) * density(17) 
  pd(39,17) = pd(39,17) + rrt(140) * density(01) 
  pd(01,01) = pd(01,01) - rrt(141) * density(18) 
  pd(01,18) = pd(01,18) - rrt(141) * density(01) 
  pd(18,01) = pd(18,01) - rrt(141) * density(18) 
  pd(18,18) = pd(18,18) - rrt(141) * density(01) 
  pd(25,01) = pd(25,01) + rrt(141) * density(18) 
  pd(25,18) = pd(25,18) + rrt(141) * density(01) 
  pd(39,01) = pd(39,01) + rrt(141) * density(18) 
  pd(39,18) = pd(39,18) + rrt(141) * density(01) 
  pd(15,01) = pd(15,01) + rrt(142) * density(16) 
  pd(15,16) = pd(15,16) + rrt(142) * density(01) 
  pd(16,01) = pd(16,01) - rrt(142) * density(16) 
  pd(16,16) = pd(16,16) - rrt(142) * density(01) 
  pd(15,01) = pd(15,01) + rrt(143) * density(17) 
  pd(15,17) = pd(15,17) + rrt(143) * density(01) 
  pd(17,01) = pd(17,01) - rrt(143) * density(17) 
  pd(17,17) = pd(17,17) - rrt(143) * density(01) 
  pd(15,01) = pd(15,01) + rrt(144) * density(18) 
  pd(15,18) = pd(15,18) + rrt(144) * density(01) 
  pd(18,01) = pd(18,01) - rrt(144) * density(18) 
  pd(18,18) = pd(18,18) - rrt(144) * density(01) 
  pd(19,01) = pd(19,01) - rrt(145) * density(19) 
  pd(19,19) = pd(19,19) - rrt(145) * density(01) 
  pd(25,01) = pd(25,01) + rrt(145) * density(19) * 2.0d0
  pd(25,19) = pd(25,19) + rrt(145) * density(01) * 2.0d0
  pd(20,01) = pd(20,01) - rrt(146) * density(20) 
  pd(20,20) = pd(20,20) - rrt(146) * density(01) 
  pd(25,01) = pd(25,01) + rrt(146) * density(20) * 2.0d0
  pd(25,20) = pd(25,20) + rrt(146) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(147) * density(19) 
  pd(01,19) = pd(01,19) + rrt(147) * density(01) 
  pd(19,01) = pd(19,01) - rrt(147) * density(19) 
  pd(19,19) = pd(19,19) - rrt(147) * density(01) 
  pd(36,01) = pd(36,01) + rrt(147) * density(19) 
  pd(36,19) = pd(36,19) + rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(20) 
  pd(01,20) = pd(01,20) + rrt(148) * density(01) 
  pd(20,01) = pd(20,01) - rrt(148) * density(20) 
  pd(20,20) = pd(20,20) - rrt(148) * density(01) 
  pd(36,01) = pd(36,01) + rrt(148) * density(20) 
  pd(36,20) = pd(36,20) + rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(19) 
  pd(01,19) = pd(01,19) + rrt(149) * density(01) 
  pd(19,01) = pd(19,01) - rrt(149) * density(19) 
  pd(19,19) = pd(19,19) - rrt(149) * density(01) 
  pd(25,01) = pd(25,01) + rrt(149) * density(19) 
  pd(25,19) = pd(25,19) + rrt(149) * density(01) 
  pd(37,01) = pd(37,01) + rrt(149) * density(19) 
  pd(37,19) = pd(37,19) + rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(20) 
  pd(01,20) = pd(01,20) + rrt(150) * density(01) 
  pd(20,01) = pd(20,01) - rrt(150) * density(20) 
  pd(20,20) = pd(20,20) - rrt(150) * density(01) 
  pd(25,01) = pd(25,01) + rrt(150) * density(20) 
  pd(25,20) = pd(25,20) + rrt(150) * density(01) 
  pd(37,01) = pd(37,01) + rrt(150) * density(20) 
  pd(37,20) = pd(37,20) + rrt(150) * density(01) 
  pd(15,01) = pd(15,01) + rrt(151) * density(24) 
  pd(15,24) = pd(15,24) + rrt(151) * density(01) 
  pd(24,01) = pd(24,01) - rrt(151) * density(24) 
  pd(24,24) = pd(24,24) - rrt(151) * density(01) 
  pd(25,01) = pd(25,01) + rrt(151) * density(24) 
  pd(25,24) = pd(25,24) + rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(24) 
  pd(01,24) = pd(01,24) + rrt(152) * density(01) 
  pd(24,01) = pd(24,01) - rrt(152) * density(24) 
  pd(24,24) = pd(24,24) - rrt(152) * density(01) 
  pd(25,01) = pd(25,01) + rrt(152) * density(24) 
  pd(25,24) = pd(25,24) + rrt(152) * density(01) 
  pd(36,01) = pd(36,01) + rrt(152) * density(24) 
  pd(36,24) = pd(36,24) + rrt(152) * density(01) 
  pd(24,01) = pd(24,01) - rrt(153) * density(24) 
  pd(24,24) = pd(24,24) - rrt(153) * density(01) 
  pd(25,01) = pd(25,01) + rrt(153) * density(24) 
  pd(25,24) = pd(25,24) + rrt(153) * density(01) 
  pd(37,01) = pd(37,01) + rrt(153) * density(24) 
  pd(37,24) = pd(37,24) + rrt(153) * density(01) 
  pd(39,01) = pd(39,01) + rrt(153) * density(24) 
  pd(39,24) = pd(39,24) + rrt(153) * density(01) 
  pd(01,01) = pd(01,01) - rrt(154) * density(24) 
  pd(01,24) = pd(01,24) - rrt(154) * density(01) 
  pd(15,01) = pd(15,01) + rrt(154) * density(24) 
  pd(15,24) = pd(15,24) + rrt(154) * density(01) 
  pd(24,01) = pd(24,01) - rrt(154) * density(24) 
  pd(24,24) = pd(24,24) - rrt(154) * density(01) 
  pd(39,01) = pd(39,01) + rrt(154) * density(24) 
  pd(39,24) = pd(39,24) + rrt(154) * density(01) 
  pd(01,01) = pd(01,01) - rrt(155) * density(24) 
  pd(01,24) = pd(01,24) - rrt(155) * density(01) 
  pd(24,01) = pd(24,01) - rrt(155) * density(24) 
  pd(24,24) = pd(24,24) - rrt(155) * density(01) 
  pd(25,01) = pd(25,01) + rrt(155) * density(24) 
  pd(25,24) = pd(25,24) + rrt(155) * density(01) 
  pd(40,01) = pd(40,01) + rrt(155) * density(24) 
  pd(40,24) = pd(40,24) + rrt(155) * density(01) 
  pd(01,01) = pd(01,01) + rrt(156) * density(25) 
  pd(01,25) = pd(01,25) + rrt(156) * density(01) 
  pd(25,01) = pd(25,01) - rrt(156) * density(25) 
  pd(25,25) = pd(25,25) - rrt(156) * density(01) 
  pd(37,01) = pd(37,01) + rrt(156) * density(25) 
  pd(37,25) = pd(37,25) + rrt(156) * density(01) 
  pd(01,01) = pd(01,01) - rrt(157) * density(26) 
  pd(01,26) = pd(01,26) - rrt(157) * density(01) 
  pd(09,01) = pd(09,01) + rrt(157) * density(26) 
  pd(09,26) = pd(09,26) + rrt(157) * density(01) 
  pd(25,01) = pd(25,01) + rrt(157) * density(26) 
  pd(25,26) = pd(25,26) + rrt(157) * density(01) 
  pd(26,01) = pd(26,01) - rrt(157) * density(26) 
  pd(26,26) = pd(26,26) - rrt(157) * density(01) 
  pd(01,01) = pd(01,01) - rrt(158) * density(26) 
  pd(01,26) = pd(01,26) - rrt(158) * density(01) 
  pd(15,01) = pd(15,01) + rrt(158) * density(26) 
  pd(15,26) = pd(15,26) + rrt(158) * density(01) 
  pd(22,01) = pd(22,01) + rrt(158) * density(26) 
  pd(22,26) = pd(22,26) + rrt(158) * density(01) 
  pd(26,01) = pd(26,01) - rrt(158) * density(26) 
  pd(26,26) = pd(26,26) - rrt(158) * density(01) 
  pd(01,01) = pd(01,01) - rrt(159) * density(27) 
  pd(01,27) = pd(01,27) - rrt(159) * density(01) 
  pd(02,01) = pd(02,01) + rrt(159) * density(27) 
  pd(02,27) = pd(02,27) + rrt(159) * density(01) 
  pd(15,01) = pd(15,01) + rrt(159) * density(27) 
  pd(15,27) = pd(15,27) + rrt(159) * density(01) 
  pd(27,01) = pd(27,01) - rrt(159) * density(27) 
  pd(27,27) = pd(27,27) - rrt(159) * density(01) 
  pd(01,01) = pd(01,01) - rrt(160) * density(28) 
  pd(01,28) = pd(01,28) - rrt(160) * density(01) 
  pd(22,01) = pd(22,01) + rrt(160) * density(28) 
  pd(22,28) = pd(22,28) + rrt(160) * density(01) 
  pd(25,01) = pd(25,01) + rrt(160) * density(28) 
  pd(25,28) = pd(25,28) + rrt(160) * density(01) 
  pd(28,01) = pd(28,01) - rrt(160) * density(28) 
  pd(28,28) = pd(28,28) - rrt(160) * density(01) 
  pd(01,01) = pd(01,01) - rrt(161) * density(29) 
  pd(01,29) = pd(01,29) - rrt(161) * density(01) 
  pd(09,01) = pd(09,01) + rrt(161) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(161) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(161) * density(29) 
  pd(29,29) = pd(29,29) - rrt(161) * density(01) 
  pd(01,01) = pd(01,01) - rrt(162) * density(30) 
  pd(01,30) = pd(01,30) - rrt(162) * density(01) 
  pd(02,01) = pd(02,01) + rrt(162) * density(30) 
  pd(02,30) = pd(02,30) + rrt(162) * density(01) 
  pd(09,01) = pd(09,01) + rrt(162) * density(30) 
  pd(09,30) = pd(09,30) + rrt(162) * density(01) 
  pd(30,01) = pd(30,01) - rrt(162) * density(30) 
  pd(30,30) = pd(30,30) - rrt(162) * density(01) 
  pd(01,01) = pd(01,01) - rrt(163) * density(31) 
  pd(01,31) = pd(01,31) - rrt(163) * density(01) 
  pd(02,01) = pd(02,01) + rrt(163) * density(31) * 2.0d0
  pd(02,31) = pd(02,31) + rrt(163) * density(01) * 2.0d0
  pd(31,01) = pd(31,01) - rrt(163) * density(31) 
  pd(31,31) = pd(31,31) - rrt(163) * density(01) 
  pd(01,01) = pd(01,01) - rrt(164) * density(33) 
  pd(01,33) = pd(01,33) - rrt(164) * density(01) 
  pd(22,01) = pd(22,01) + rrt(164) * density(33) * 2.0d0
  pd(22,33) = pd(22,33) + rrt(164) * density(01) * 2.0d0
  pd(33,01) = pd(33,01) - rrt(164) * density(33) 
  pd(33,33) = pd(33,33) - rrt(164) * density(01) 
  pd(01,01) = pd(01,01) - rrt(165) * density(15) 
  pd(01,15) = pd(01,15) - rrt(165) * density(01) 
  pd(15,01) = pd(15,01) - rrt(165) * density(15) 
  pd(15,15) = pd(15,15) - rrt(165) * density(01) 
  pd(40,01) = pd(40,01) + rrt(165) * density(15) 
  pd(40,15) = pd(40,15) + rrt(165) * density(01) 
  pd(01,01) = pd(01,01) - rrt(166) * density(24) 
  pd(01,24) = pd(01,24) - rrt(166) * density(01) 
  pd(24,01) = pd(24,01) - rrt(166) * density(24) 
  pd(24,24) = pd(24,24) - rrt(166) * density(01) 
  pd(41,01) = pd(41,01) + rrt(166) * density(24) 
  pd(41,24) = pd(41,24) + rrt(166) * density(01) 
  pd(01,01) = pd(01,01) - rrt(167) * density(25) 
  pd(01,25) = pd(01,25) - rrt(167) * density(01) 
  pd(25,01) = pd(25,01) - rrt(167) * density(25) 
  pd(25,25) = pd(25,25) - rrt(167) * density(01) 
  pd(39,01) = pd(39,01) + rrt(167) * density(25) 
  pd(39,25) = pd(39,25) + rrt(167) * density(01) 
  pd(01,01) = pd(01,01) - rrt(168) * density(36) 
  pd(01,36) = pd(01,36) - rrt(168) * density(01) 
  pd(15,01) = pd(15,01) + rrt(168) * density(36) 
  pd(15,36) = pd(15,36) + rrt(168) * density(01) 
  pd(36,01) = pd(36,01) - rrt(168) * density(36) 
  pd(36,36) = pd(36,36) - rrt(168) * density(01) 
  pd(01,01) = pd(01,01) - rrt(169) * density(36) 
  pd(01,36) = pd(01,36) - rrt(169) * density(01) 
  pd(25,01) = pd(25,01) + rrt(169) * density(36) * 2.0d0
  pd(25,36) = pd(25,36) + rrt(169) * density(01) * 2.0d0
  pd(36,01) = pd(36,01) - rrt(169) * density(36) 
  pd(36,36) = pd(36,36) - rrt(169) * density(01) 
  pd(01,01) = pd(01,01) - rrt(170) * density(37) 
  pd(01,37) = pd(01,37) - rrt(170) * density(01) 
  pd(25,01) = pd(25,01) + rrt(170) * density(37) 
  pd(25,37) = pd(25,37) + rrt(170) * density(01) 
  pd(37,01) = pd(37,01) - rrt(170) * density(37) 
  pd(37,37) = pd(37,37) - rrt(170) * density(01) 
  pd(01,01) = pd(01,01) - rrt(171) * density(38) 
  pd(01,38) = pd(01,38) - rrt(171) * density(01) 
  pd(15,01) = pd(15,01) + rrt(171) * density(38) * 2.0d0
  pd(15,38) = pd(15,38) + rrt(171) * density(01) * 2.0d0
  pd(38,01) = pd(38,01) - rrt(171) * density(38) 
  pd(38,38) = pd(38,38) - rrt(171) * density(01) 
  pd(02,02) = pd(02,02) - rrt(172) * density(25) 
  pd(02,25) = pd(02,25) - rrt(172) * density(02) 
  pd(09,02) = pd(09,02) + rrt(172) * density(25) 
  pd(09,25) = pd(09,25) + rrt(172) * density(02) 
  pd(15,02) = pd(15,02) + rrt(172) * density(25) 
  pd(15,25) = pd(15,25) + rrt(172) * density(02) 
  pd(25,02) = pd(25,02) - rrt(172) * density(25) 
  pd(25,25) = pd(25,25) - rrt(172) * density(02) 
  pd(02,02) = pd(02,02) - rrt(173) * density(22) 
  pd(02,22) = pd(02,22) - rrt(173) * density(02) 
  pd(09,02) = pd(09,02) + rrt(173) * density(22) * 2.0d0
  pd(09,22) = pd(09,22) + rrt(173) * density(02) * 2.0d0
  pd(22,02) = pd(22,02) - rrt(173) * density(22) 
  pd(22,22) = pd(22,22) - rrt(173) * density(02) 
  pd(03,03) = pd(03,03) - rrt(174) * density(25) 
  pd(03,25) = pd(03,25) - rrt(174) * density(03) 
  pd(09,03) = pd(09,03) + rrt(174) * density(25) 
  pd(09,25) = pd(09,25) + rrt(174) * density(03) 
  pd(15,03) = pd(15,03) + rrt(174) * density(25) 
  pd(15,25) = pd(15,25) + rrt(174) * density(03) 
  pd(25,03) = pd(25,03) - rrt(174) * density(25) 
  pd(25,25) = pd(25,25) - rrt(174) * density(03) 
  pd(04,04) = pd(04,04) - rrt(175) * density(25) 
  pd(04,25) = pd(04,25) - rrt(175) * density(04) 
  pd(09,04) = pd(09,04) + rrt(175) * density(25) 
  pd(09,25) = pd(09,25) + rrt(175) * density(04) 
  pd(15,04) = pd(15,04) + rrt(175) * density(25) 
  pd(15,25) = pd(15,25) + rrt(175) * density(04) 
  pd(25,04) = pd(25,04) - rrt(175) * density(25) 
  pd(25,25) = pd(25,25) - rrt(175) * density(04) 
  pd(05,05) = pd(05,05) - rrt(176) * density(25) 
  pd(05,25) = pd(05,25) - rrt(176) * density(05) 
  pd(09,05) = pd(09,05) + rrt(176) * density(25) 
  pd(09,25) = pd(09,25) + rrt(176) * density(05) 
  pd(15,05) = pd(15,05) + rrt(176) * density(25) 
  pd(15,25) = pd(15,25) + rrt(176) * density(05) 
  pd(25,05) = pd(25,05) - rrt(176) * density(25) 
  pd(25,25) = pd(25,25) - rrt(176) * density(05) 
  pd(06,06) = pd(06,06) - rrt(177) * density(25) 
  pd(06,25) = pd(06,25) - rrt(177) * density(06) 
  pd(09,06) = pd(09,06) + rrt(177) * density(25) 
  pd(09,25) = pd(09,25) + rrt(177) * density(06) 
  pd(15,06) = pd(15,06) + rrt(177) * density(25) 
  pd(15,25) = pd(15,25) + rrt(177) * density(06) 
  pd(25,06) = pd(25,06) - rrt(177) * density(25) 
  pd(25,25) = pd(25,25) - rrt(177) * density(06) 
  pd(03,03) = pd(03,03) - rrt(178) * density(22) 
  pd(03,22) = pd(03,22) - rrt(178) * density(03) 
  pd(09,03) = pd(09,03) + rrt(178) * density(22) * 2.0d0
  pd(09,22) = pd(09,22) + rrt(178) * density(03) * 2.0d0
  pd(22,03) = pd(22,03) - rrt(178) * density(22) 
  pd(22,22) = pd(22,22) - rrt(178) * density(03) 
  pd(04,04) = pd(04,04) - rrt(179) * density(22) 
  pd(04,22) = pd(04,22) - rrt(179) * density(04) 
  pd(09,04) = pd(09,04) + rrt(179) * density(22) * 2.0d0
  pd(09,22) = pd(09,22) + rrt(179) * density(04) * 2.0d0
  pd(22,04) = pd(22,04) - rrt(179) * density(22) 
  pd(22,22) = pd(22,22) - rrt(179) * density(04) 
  pd(05,05) = pd(05,05) - rrt(180) * density(22) 
  pd(05,22) = pd(05,22) - rrt(180) * density(05) 
  pd(09,05) = pd(09,05) + rrt(180) * density(22) * 2.0d0
  pd(09,22) = pd(09,22) + rrt(180) * density(05) * 2.0d0
  pd(22,05) = pd(22,05) - rrt(180) * density(22) 
  pd(22,22) = pd(22,22) - rrt(180) * density(05) 
  pd(06,06) = pd(06,06) - rrt(181) * density(22) 
  pd(06,22) = pd(06,22) - rrt(181) * density(06) 
  pd(09,06) = pd(09,06) + rrt(181) * density(22) * 2.0d0
  pd(09,22) = pd(09,22) + rrt(181) * density(06) * 2.0d0
  pd(22,06) = pd(22,06) - rrt(181) * density(22) 
  pd(22,22) = pd(22,22) - rrt(181) * density(06) 
  pd(02,09) = pd(02,09) + rrt(182) * density(25) 
  pd(02,25) = pd(02,25) + rrt(182) * density(09) 
  pd(09,09) = pd(09,09) - rrt(182) * density(25) 
  pd(09,25) = pd(09,25) - rrt(182) * density(09) 
  pd(25,09) = pd(25,09) - rrt(182) * density(25) 
  pd(25,25) = pd(25,25) - rrt(182) * density(09) 
  pd(02,09) = pd(02,09) + rrt(183) * density(24) 
  pd(02,24) = pd(02,24) + rrt(183) * density(09) 
  pd(09,09) = pd(09,09) - rrt(183) * density(24) 
  pd(09,24) = pd(09,24) - rrt(183) * density(09) 
  pd(15,09) = pd(15,09) + rrt(183) * density(24) 
  pd(15,24) = pd(15,24) + rrt(183) * density(09) 
  pd(24,09) = pd(24,09) - rrt(183) * density(24) 
  pd(24,24) = pd(24,24) - rrt(183) * density(09) 
  pd(09,09) = pd(09,09) - rrt(184) * density(22) 
  pd(09,22) = pd(09,22) - rrt(184) * density(09) 
  pd(21,09) = pd(21,09) + rrt(184) * density(22) 
  pd(21,22) = pd(21,22) + rrt(184) * density(09) 
  pd(22,09) = pd(22,09) - rrt(184) * density(22) 
  pd(22,22) = pd(22,22) - rrt(184) * density(09) 
  pd(02,09) = pd(02,09) + rrt(185) * density(15) 
  pd(02,15) = pd(02,15) + rrt(185) * density(09) 
  pd(09,09) = pd(09,09) - rrt(185) * density(15) 
  pd(09,15) = pd(09,15) - rrt(185) * density(09) 
  pd(15,09) = pd(15,09) - rrt(185) * density(15) 
  pd(15,15) = pd(15,15) - rrt(185) * density(09) 
  pd(25,09) = pd(25,09) + rrt(185) * density(15) 
  pd(25,15) = pd(25,15) + rrt(185) * density(09) 
  pd(02,10) = pd(02,10) + rrt(186) * density(25) 
  pd(02,25) = pd(02,25) + rrt(186) * density(10) 
  pd(10,10) = pd(10,10) - rrt(186) * density(25) 
  pd(10,25) = pd(10,25) - rrt(186) * density(10) 
  pd(25,10) = pd(25,10) - rrt(186) * density(25) 
  pd(25,25) = pd(25,25) - rrt(186) * density(10) 
  pd(02,09) = pd(02,09) + rrt(187) * density(16) 
  pd(02,16) = pd(02,16) + rrt(187) * density(09) 
  pd(09,09) = pd(09,09) - rrt(187) * density(16) 
  pd(09,16) = pd(09,16) - rrt(187) * density(09) 
  pd(16,09) = pd(16,09) - rrt(187) * density(16) 
  pd(16,16) = pd(16,16) - rrt(187) * density(09) 
  pd(25,09) = pd(25,09) + rrt(187) * density(16) 
  pd(25,16) = pd(25,16) + rrt(187) * density(09) 
  pd(02,09) = pd(02,09) + rrt(188) * density(17) 
  pd(02,17) = pd(02,17) + rrt(188) * density(09) 
  pd(09,09) = pd(09,09) - rrt(188) * density(17) 
  pd(09,17) = pd(09,17) - rrt(188) * density(09) 
  pd(17,09) = pd(17,09) - rrt(188) * density(17) 
  pd(17,17) = pd(17,17) - rrt(188) * density(09) 
  pd(25,09) = pd(25,09) + rrt(188) * density(17) 
  pd(25,17) = pd(25,17) + rrt(188) * density(09) 
  pd(02,09) = pd(02,09) + rrt(189) * density(18) 
  pd(02,18) = pd(02,18) + rrt(189) * density(09) 
  pd(09,09) = pd(09,09) - rrt(189) * density(18) 
  pd(09,18) = pd(09,18) - rrt(189) * density(09) 
  pd(18,09) = pd(18,09) - rrt(189) * density(18) 
  pd(18,18) = pd(18,18) - rrt(189) * density(09) 
  pd(25,09) = pd(25,09) + rrt(189) * density(18) 
  pd(25,18) = pd(25,18) + rrt(189) * density(09) 
  pd(02,10) = pd(02,10) + rrt(190) * density(15) 
  pd(02,15) = pd(02,15) + rrt(190) * density(10) 
  pd(10,10) = pd(10,10) - rrt(190) * density(15) 
  pd(10,15) = pd(10,15) - rrt(190) * density(10) 
  pd(15,10) = pd(15,10) - rrt(190) * density(15) 
  pd(15,15) = pd(15,15) - rrt(190) * density(10) 
  pd(25,10) = pd(25,10) + rrt(190) * density(15) 
  pd(25,15) = pd(25,15) + rrt(190) * density(10) 
  pd(02,10) = pd(02,10) + rrt(191) * density(16) 
  pd(02,16) = pd(02,16) + rrt(191) * density(10) 
  pd(10,10) = pd(10,10) - rrt(191) * density(16) 
  pd(10,16) = pd(10,16) - rrt(191) * density(10) 
  pd(16,10) = pd(16,10) - rrt(191) * density(16) 
  pd(16,16) = pd(16,16) - rrt(191) * density(10) 
  pd(25,10) = pd(25,10) + rrt(191) * density(16) 
  pd(25,16) = pd(25,16) + rrt(191) * density(10) 
  pd(02,10) = pd(02,10) + rrt(192) * density(17) 
  pd(02,17) = pd(02,17) + rrt(192) * density(10) 
  pd(10,10) = pd(10,10) - rrt(192) * density(17) 
  pd(10,17) = pd(10,17) - rrt(192) * density(10) 
  pd(17,10) = pd(17,10) - rrt(192) * density(17) 
  pd(17,17) = pd(17,17) - rrt(192) * density(10) 
  pd(25,10) = pd(25,10) + rrt(192) * density(17) 
  pd(25,17) = pd(25,17) + rrt(192) * density(10) 
  pd(02,10) = pd(02,10) + rrt(193) * density(18) 
  pd(02,18) = pd(02,18) + rrt(193) * density(10) 
  pd(10,10) = pd(10,10) - rrt(193) * density(18) 
  pd(10,18) = pd(10,18) - rrt(193) * density(10) 
  pd(18,10) = pd(18,10) - rrt(193) * density(18) 
  pd(18,18) = pd(18,18) - rrt(193) * density(10) 
  pd(25,10) = pd(25,10) + rrt(193) * density(18) 
  pd(25,18) = pd(25,18) + rrt(193) * density(10) 
  pd(02,10) = pd(02,10) + rrt(194) * density(24) 
  pd(02,24) = pd(02,24) + rrt(194) * density(10) 
  pd(10,10) = pd(10,10) - rrt(194) * density(24) 
  pd(10,24) = pd(10,24) - rrt(194) * density(10) 
  pd(15,10) = pd(15,10) + rrt(194) * density(24) 
  pd(15,24) = pd(15,24) + rrt(194) * density(10) 
  pd(24,10) = pd(24,10) - rrt(194) * density(24) 
  pd(24,24) = pd(24,24) - rrt(194) * density(10) 
  pd(10,10) = pd(10,10) - rrt(195) * density(22) 
  pd(10,22) = pd(10,22) - rrt(195) * density(10) 
  pd(21,10) = pd(21,10) + rrt(195) * density(22) 
  pd(21,22) = pd(21,22) + rrt(195) * density(10) 
  pd(22,10) = pd(22,10) - rrt(195) * density(22) 
  pd(22,22) = pd(22,22) - rrt(195) * density(10) 
  pd(09,15) = pd(09,15) + rrt(196) * density(22) 
  pd(09,22) = pd(09,22) + rrt(196) * density(15) 
  pd(15,15) = pd(15,15) - rrt(196) * density(22) 
  pd(15,22) = pd(15,22) - rrt(196) * density(15) 
  pd(22,15) = pd(22,15) - rrt(196) * density(22) 
  pd(22,22) = pd(22,22) - rrt(196) * density(15) 
  pd(25,15) = pd(25,15) + rrt(196) * density(22) 
  pd(25,22) = pd(25,22) + rrt(196) * density(15) 
  pd(02,15) = pd(02,15) + rrt(197) * density(21) 
  pd(02,21) = pd(02,21) + rrt(197) * density(15) 
  pd(09,15) = pd(09,15) + rrt(197) * density(21) 
  pd(09,21) = pd(09,21) + rrt(197) * density(15) 
  pd(15,15) = pd(15,15) - rrt(197) * density(21) 
  pd(15,21) = pd(15,21) - rrt(197) * density(15) 
  pd(21,15) = pd(21,15) - rrt(197) * density(21) 
  pd(21,21) = pd(21,21) - rrt(197) * density(15) 
  pd(15,15) = pd(15,15) - rrt(198) * density(25) 
  pd(15,25) = pd(15,25) - rrt(198) * density(15) 
  pd(24,15) = pd(24,15) + rrt(198) * density(25) 
  pd(24,25) = pd(24,25) + rrt(198) * density(15) 
  pd(25,15) = pd(25,15) - rrt(198) * density(25) 
  pd(25,25) = pd(25,25) - rrt(198) * density(15) 
  pd(09,16) = pd(09,16) + rrt(199) * density(22) 
  pd(09,22) = pd(09,22) + rrt(199) * density(16) 
  pd(16,16) = pd(16,16) - rrt(199) * density(22) 
  pd(16,22) = pd(16,22) - rrt(199) * density(16) 
  pd(22,16) = pd(22,16) - rrt(199) * density(22) 
  pd(22,22) = pd(22,22) - rrt(199) * density(16) 
  pd(25,16) = pd(25,16) + rrt(199) * density(22) 
  pd(25,22) = pd(25,22) + rrt(199) * density(16) 
  pd(09,17) = pd(09,17) + rrt(200) * density(22) 
  pd(09,22) = pd(09,22) + rrt(200) * density(17) 
  pd(17,17) = pd(17,17) - rrt(200) * density(22) 
  pd(17,22) = pd(17,22) - rrt(200) * density(17) 
  pd(22,17) = pd(22,17) - rrt(200) * density(22) 
  pd(22,22) = pd(22,22) - rrt(200) * density(17) 
  pd(25,17) = pd(25,17) + rrt(200) * density(22) 
  pd(25,22) = pd(25,22) + rrt(200) * density(17) 
  pd(09,18) = pd(09,18) + rrt(201) * density(22) 
  pd(09,22) = pd(09,22) + rrt(201) * density(18) 
  pd(18,18) = pd(18,18) - rrt(201) * density(22) 
  pd(18,22) = pd(18,22) - rrt(201) * density(18) 
  pd(22,18) = pd(22,18) - rrt(201) * density(22) 
  pd(22,22) = pd(22,22) - rrt(201) * density(18) 
  pd(25,18) = pd(25,18) + rrt(201) * density(22) 
  pd(25,22) = pd(25,22) + rrt(201) * density(18) 
  pd(16,16) = pd(16,16) - rrt(202) * density(25) 
  pd(16,25) = pd(16,25) - rrt(202) * density(16) 
  pd(24,16) = pd(24,16) + rrt(202) * density(25) 
  pd(24,25) = pd(24,25) + rrt(202) * density(16) 
  pd(25,16) = pd(25,16) - rrt(202) * density(25) 
  pd(25,25) = pd(25,25) - rrt(202) * density(16) 
  pd(17,17) = pd(17,17) - rrt(203) * density(25) 
  pd(17,25) = pd(17,25) - rrt(203) * density(17) 
  pd(24,17) = pd(24,17) + rrt(203) * density(25) 
  pd(24,25) = pd(24,25) + rrt(203) * density(17) 
  pd(25,17) = pd(25,17) - rrt(203) * density(25) 
  pd(25,25) = pd(25,25) - rrt(203) * density(17) 
  pd(18,18) = pd(18,18) - rrt(204) * density(25) 
  pd(18,25) = pd(18,25) - rrt(204) * density(18) 
  pd(24,18) = pd(24,18) + rrt(204) * density(25) 
  pd(24,25) = pd(24,25) + rrt(204) * density(18) 
  pd(25,18) = pd(25,18) - rrt(204) * density(25) 
  pd(25,25) = pd(25,25) - rrt(204) * density(18) 
  pd(09,22) = pd(09,22) + rrt(205) * density(25) 
  pd(09,25) = pd(09,25) + rrt(205) * density(22) 
  pd(22,22) = pd(22,22) - rrt(205) * density(25) 
  pd(22,25) = pd(22,25) - rrt(205) * density(22) 
  pd(25,22) = pd(25,22) - rrt(205) * density(25) 
  pd(25,25) = pd(25,25) - rrt(205) * density(22) 
  pd(09,21) = pd(09,21) + rrt(206) * density(25) * 2.0d0
  pd(09,25) = pd(09,25) + rrt(206) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(206) * density(25) 
  pd(21,25) = pd(21,25) - rrt(206) * density(21) 
  pd(25,21) = pd(25,21) - rrt(206) * density(25) 
  pd(25,25) = pd(25,25) - rrt(206) * density(21) 
  pd(15,24) = pd(15,24) + rrt(207) * density(25) * 2.0d0
  pd(15,25) = pd(15,25) + rrt(207) * density(24) * 2.0d0
  pd(24,24) = pd(24,24) - rrt(207) * density(25) 
  pd(24,25) = pd(24,25) - rrt(207) * density(24) 
  pd(25,24) = pd(25,24) - rrt(207) * density(25) 
  pd(25,25) = pd(25,25) - rrt(207) * density(24) 
  pd(15,24) = pd(15,24) + rrt(208) 
  pd(24,24) = pd(24,24) - rrt(208) 
  pd(25,24) = pd(25,24) + rrt(208) 
  pd(15,25) = pd(15,25) + rrt(209) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(209) * density(25) * 4.0d0
  pd(02,02) = pd(02,02) - rrt(210) * density(36) 
  pd(02,36) = pd(02,36) - rrt(210) * density(02) 
  pd(27,02) = pd(27,02) + rrt(210) * density(36) 
  pd(27,36) = pd(27,36) + rrt(210) * density(02) 
  pd(36,02) = pd(36,02) - rrt(210) * density(36) 
  pd(36,36) = pd(36,36) - rrt(210) * density(02) 
  pd(02,02) = pd(02,02) - rrt(211) * density(37) 
  pd(02,37) = pd(02,37) - rrt(211) * density(02) 
  pd(09,02) = pd(09,02) + rrt(211) * density(37) 
  pd(09,37) = pd(09,37) + rrt(211) * density(02) 
  pd(36,02) = pd(36,02) + rrt(211) * density(37) 
  pd(36,37) = pd(36,37) + rrt(211) * density(02) 
  pd(37,02) = pd(37,02) - rrt(211) * density(37) 
  pd(37,37) = pd(37,37) - rrt(211) * density(02) 
  pd(02,02) = pd(02,02) - rrt(212) * density(37) 
  pd(02,37) = pd(02,37) - rrt(212) * density(02) 
  pd(25,02) = pd(25,02) + rrt(212) * density(37) 
  pd(25,37) = pd(25,37) + rrt(212) * density(02) 
  pd(26,02) = pd(26,02) + rrt(212) * density(37) 
  pd(26,37) = pd(26,37) + rrt(212) * density(02) 
  pd(37,02) = pd(37,02) - rrt(212) * density(37) 
  pd(37,37) = pd(37,37) - rrt(212) * density(02) 
  pd(02,02) = pd(02,02) - rrt(213) * density(32) 
  pd(02,32) = pd(02,32) - rrt(213) * density(02) 
  pd(09,02) = pd(09,02) + rrt(213) * density(32) 
  pd(09,32) = pd(09,32) + rrt(213) * density(02) 
  pd(28,02) = pd(28,02) + rrt(213) * density(32) 
  pd(28,32) = pd(28,32) + rrt(213) * density(02) 
  pd(32,02) = pd(32,02) - rrt(213) * density(32) 
  pd(32,32) = pd(32,32) - rrt(213) * density(02) 
  pd(02,02) = pd(02,02) - rrt(214) * density(28) 
  pd(02,28) = pd(02,28) - rrt(214) * density(02) 
  pd(09,02) = pd(09,02) + rrt(214) * density(28) 
  pd(09,28) = pd(09,28) + rrt(214) * density(02) 
  pd(26,02) = pd(26,02) + rrt(214) * density(28) 
  pd(26,28) = pd(26,28) + rrt(214) * density(02) 
  pd(28,02) = pd(28,02) - rrt(214) * density(28) 
  pd(28,28) = pd(28,28) - rrt(214) * density(02) 
  pd(02,02) = pd(02,02) - rrt(215) * density(39) 
  pd(02,39) = pd(02,39) - rrt(215) * density(02) 
  pd(34,02) = pd(34,02) + rrt(215) * density(39) 
  pd(34,39) = pd(34,39) + rrt(215) * density(02) 
  pd(39,02) = pd(39,02) - rrt(215) * density(39) 
  pd(39,39) = pd(39,39) - rrt(215) * density(02) 
  pd(02,02) = pd(02,02) - rrt(216) * density(40) 
  pd(02,40) = pd(02,40) - rrt(216) * density(02) 
  pd(35,02) = pd(35,02) + rrt(216) * density(40) 
  pd(35,40) = pd(35,40) + rrt(216) * density(02) 
  pd(40,02) = pd(40,02) - rrt(216) * density(40) 
  pd(40,40) = pd(40,40) - rrt(216) * density(02) 
  pd(02,02) = pd(02,02) - rrt(217) * density(41) 
  pd(02,41) = pd(02,41) - rrt(217) * density(02) 
  pd(15,02) = pd(15,02) + rrt(217) * density(41) 
  pd(15,41) = pd(15,41) + rrt(217) * density(02) 
  pd(34,02) = pd(34,02) + rrt(217) * density(41) 
  pd(34,41) = pd(34,41) + rrt(217) * density(02) 
  pd(41,02) = pd(41,02) - rrt(217) * density(41) 
  pd(41,41) = pd(41,41) - rrt(217) * density(02) 
  pd(02,02) = pd(02,02) - rrt(218) * density(42) 
  pd(02,42) = pd(02,42) - rrt(218) * density(02) 
  pd(15,02) = pd(15,02) + rrt(218) * density(42) 
  pd(15,42) = pd(15,42) + rrt(218) * density(02) 
  pd(35,02) = pd(35,02) + rrt(218) * density(42) 
  pd(35,42) = pd(35,42) + rrt(218) * density(02) 
  pd(42,02) = pd(42,02) - rrt(218) * density(42) 
  pd(42,42) = pd(42,42) - rrt(218) * density(02) 
  pd(02,02) = pd(02,02) - rrt(219) * density(26) 
  pd(02,26) = pd(02,26) - rrt(219) * density(02) 
  pd(26,02) = pd(26,02) - rrt(219) * density(26) 
  pd(26,26) = pd(26,26) - rrt(219) * density(02) 
  pd(31,02) = pd(31,02) + rrt(219) * density(26) 
  pd(31,26) = pd(31,26) + rrt(219) * density(02) 
  pd(03,03) = pd(03,03) - rrt(220) * density(37) 
  pd(03,37) = pd(03,37) - rrt(220) * density(03) 
  pd(09,03) = pd(09,03) + rrt(220) * density(37) 
  pd(09,37) = pd(09,37) + rrt(220) * density(03) 
  pd(36,03) = pd(36,03) + rrt(220) * density(37) 
  pd(36,37) = pd(36,37) + rrt(220) * density(03) 
  pd(37,03) = pd(37,03) - rrt(220) * density(37) 
  pd(37,37) = pd(37,37) - rrt(220) * density(03) 
  pd(04,04) = pd(04,04) - rrt(221) * density(37) 
  pd(04,37) = pd(04,37) - rrt(221) * density(04) 
  pd(09,04) = pd(09,04) + rrt(221) * density(37) 
  pd(09,37) = pd(09,37) + rrt(221) * density(04) 
  pd(36,04) = pd(36,04) + rrt(221) * density(37) 
  pd(36,37) = pd(36,37) + rrt(221) * density(04) 
  pd(37,04) = pd(37,04) - rrt(221) * density(37) 
  pd(37,37) = pd(37,37) - rrt(221) * density(04) 
  pd(05,05) = pd(05,05) - rrt(222) * density(37) 
  pd(05,37) = pd(05,37) - rrt(222) * density(05) 
  pd(09,05) = pd(09,05) + rrt(222) * density(37) 
  pd(09,37) = pd(09,37) + rrt(222) * density(05) 
  pd(36,05) = pd(36,05) + rrt(222) * density(37) 
  pd(36,37) = pd(36,37) + rrt(222) * density(05) 
  pd(37,05) = pd(37,05) - rrt(222) * density(37) 
  pd(37,37) = pd(37,37) - rrt(222) * density(05) 
  pd(06,06) = pd(06,06) - rrt(223) * density(37) 
  pd(06,37) = pd(06,37) - rrt(223) * density(06) 
  pd(09,06) = pd(09,06) + rrt(223) * density(37) 
  pd(09,37) = pd(09,37) + rrt(223) * density(06) 
  pd(36,06) = pd(36,06) + rrt(223) * density(37) 
  pd(36,37) = pd(36,37) + rrt(223) * density(06) 
  pd(37,06) = pd(37,06) - rrt(223) * density(37) 
  pd(37,37) = pd(37,37) - rrt(223) * density(06) 
  pd(03,03) = pd(03,03) - rrt(224) * density(37) 
  pd(03,37) = pd(03,37) - rrt(224) * density(03) 
  pd(25,03) = pd(25,03) + rrt(224) * density(37) 
  pd(25,37) = pd(25,37) + rrt(224) * density(03) 
  pd(26,03) = pd(26,03) + rrt(224) * density(37) 
  pd(26,37) = pd(26,37) + rrt(224) * density(03) 
  pd(37,03) = pd(37,03) - rrt(224) * density(37) 
  pd(37,37) = pd(37,37) - rrt(224) * density(03) 
  pd(04,04) = pd(04,04) - rrt(225) * density(37) 
  pd(04,37) = pd(04,37) - rrt(225) * density(04) 
  pd(25,04) = pd(25,04) + rrt(225) * density(37) 
  pd(25,37) = pd(25,37) + rrt(225) * density(04) 
  pd(26,04) = pd(26,04) + rrt(225) * density(37) 
  pd(26,37) = pd(26,37) + rrt(225) * density(04) 
  pd(37,04) = pd(37,04) - rrt(225) * density(37) 
  pd(37,37) = pd(37,37) - rrt(225) * density(04) 
  pd(05,05) = pd(05,05) - rrt(226) * density(37) 
  pd(05,37) = pd(05,37) - rrt(226) * density(05) 
  pd(25,05) = pd(25,05) + rrt(226) * density(37) 
  pd(25,37) = pd(25,37) + rrt(226) * density(05) 
  pd(26,05) = pd(26,05) + rrt(226) * density(37) 
  pd(26,37) = pd(26,37) + rrt(226) * density(05) 
  pd(37,05) = pd(37,05) - rrt(226) * density(37) 
  pd(37,37) = pd(37,37) - rrt(226) * density(05) 
  pd(06,06) = pd(06,06) - rrt(227) * density(37) 
  pd(06,37) = pd(06,37) - rrt(227) * density(06) 
  pd(25,06) = pd(25,06) + rrt(227) * density(37) 
  pd(25,37) = pd(25,37) + rrt(227) * density(06) 
  pd(26,06) = pd(26,06) + rrt(227) * density(37) 
  pd(26,37) = pd(26,37) + rrt(227) * density(06) 
  pd(37,06) = pd(37,06) - rrt(227) * density(37) 
  pd(37,37) = pd(37,37) - rrt(227) * density(06) 
  pd(03,03) = pd(03,03) - rrt(228) * density(32) 
  pd(03,32) = pd(03,32) - rrt(228) * density(03) 
  pd(09,03) = pd(09,03) + rrt(228) * density(32) 
  pd(09,32) = pd(09,32) + rrt(228) * density(03) 
  pd(28,03) = pd(28,03) + rrt(228) * density(32) 
  pd(28,32) = pd(28,32) + rrt(228) * density(03) 
  pd(32,03) = pd(32,03) - rrt(228) * density(32) 
  pd(32,32) = pd(32,32) - rrt(228) * density(03) 
  pd(04,04) = pd(04,04) - rrt(229) * density(32) 
  pd(04,32) = pd(04,32) - rrt(229) * density(04) 
  pd(09,04) = pd(09,04) + rrt(229) * density(32) 
  pd(09,32) = pd(09,32) + rrt(229) * density(04) 
  pd(28,04) = pd(28,04) + rrt(229) * density(32) 
  pd(28,32) = pd(28,32) + rrt(229) * density(04) 
  pd(32,04) = pd(32,04) - rrt(229) * density(32) 
  pd(32,32) = pd(32,32) - rrt(229) * density(04) 
  pd(05,05) = pd(05,05) - rrt(230) * density(32) 
  pd(05,32) = pd(05,32) - rrt(230) * density(05) 
  pd(09,05) = pd(09,05) + rrt(230) * density(32) 
  pd(09,32) = pd(09,32) + rrt(230) * density(05) 
  pd(28,05) = pd(28,05) + rrt(230) * density(32) 
  pd(28,32) = pd(28,32) + rrt(230) * density(05) 
  pd(32,05) = pd(32,05) - rrt(230) * density(32) 
  pd(32,32) = pd(32,32) - rrt(230) * density(05) 
  pd(06,06) = pd(06,06) - rrt(231) * density(32) 
  pd(06,32) = pd(06,32) - rrt(231) * density(06) 
  pd(09,06) = pd(09,06) + rrt(231) * density(32) 
  pd(09,32) = pd(09,32) + rrt(231) * density(06) 
  pd(28,06) = pd(28,06) + rrt(231) * density(32) 
  pd(28,32) = pd(28,32) + rrt(231) * density(06) 
  pd(32,06) = pd(32,06) - rrt(231) * density(32) 
  pd(32,32) = pd(32,32) - rrt(231) * density(06) 
  pd(03,03) = pd(03,03) - rrt(232) * density(28) 
  pd(03,28) = pd(03,28) - rrt(232) * density(03) 
  pd(09,03) = pd(09,03) + rrt(232) * density(28) 
  pd(09,28) = pd(09,28) + rrt(232) * density(03) 
  pd(26,03) = pd(26,03) + rrt(232) * density(28) 
  pd(26,28) = pd(26,28) + rrt(232) * density(03) 
  pd(28,03) = pd(28,03) - rrt(232) * density(28) 
  pd(28,28) = pd(28,28) - rrt(232) * density(03) 
  pd(04,04) = pd(04,04) - rrt(233) * density(28) 
  pd(04,28) = pd(04,28) - rrt(233) * density(04) 
  pd(09,04) = pd(09,04) + rrt(233) * density(28) 
  pd(09,28) = pd(09,28) + rrt(233) * density(04) 
  pd(26,04) = pd(26,04) + rrt(233) * density(28) 
  pd(26,28) = pd(26,28) + rrt(233) * density(04) 
  pd(28,04) = pd(28,04) - rrt(233) * density(28) 
  pd(28,28) = pd(28,28) - rrt(233) * density(04) 
  pd(05,05) = pd(05,05) - rrt(234) * density(28) 
  pd(05,28) = pd(05,28) - rrt(234) * density(05) 
  pd(09,05) = pd(09,05) + rrt(234) * density(28) 
  pd(09,28) = pd(09,28) + rrt(234) * density(05) 
  pd(26,05) = pd(26,05) + rrt(234) * density(28) 
  pd(26,28) = pd(26,28) + rrt(234) * density(05) 
  pd(28,05) = pd(28,05) - rrt(234) * density(28) 
  pd(28,28) = pd(28,28) - rrt(234) * density(05) 
  pd(06,06) = pd(06,06) - rrt(235) * density(28) 
  pd(06,28) = pd(06,28) - rrt(235) * density(06) 
  pd(09,06) = pd(09,06) + rrt(235) * density(28) 
  pd(09,28) = pd(09,28) + rrt(235) * density(06) 
  pd(26,06) = pd(26,06) + rrt(235) * density(28) 
  pd(26,28) = pd(26,28) + rrt(235) * density(06) 
  pd(28,06) = pd(28,06) - rrt(235) * density(28) 
  pd(28,28) = pd(28,28) - rrt(235) * density(06) 
  pd(03,03) = pd(03,03) - rrt(236) * density(39) 
  pd(03,39) = pd(03,39) - rrt(236) * density(03) 
  pd(34,03) = pd(34,03) + rrt(236) * density(39) 
  pd(34,39) = pd(34,39) + rrt(236) * density(03) 
  pd(39,03) = pd(39,03) - rrt(236) * density(39) 
  pd(39,39) = pd(39,39) - rrt(236) * density(03) 
  pd(04,04) = pd(04,04) - rrt(237) * density(39) 
  pd(04,39) = pd(04,39) - rrt(237) * density(04) 
  pd(34,04) = pd(34,04) + rrt(237) * density(39) 
  pd(34,39) = pd(34,39) + rrt(237) * density(04) 
  pd(39,04) = pd(39,04) - rrt(237) * density(39) 
  pd(39,39) = pd(39,39) - rrt(237) * density(04) 
  pd(05,05) = pd(05,05) - rrt(238) * density(39) 
  pd(05,39) = pd(05,39) - rrt(238) * density(05) 
  pd(34,05) = pd(34,05) + rrt(238) * density(39) 
  pd(34,39) = pd(34,39) + rrt(238) * density(05) 
  pd(39,05) = pd(39,05) - rrt(238) * density(39) 
  pd(39,39) = pd(39,39) - rrt(238) * density(05) 
  pd(06,06) = pd(06,06) - rrt(239) * density(39) 
  pd(06,39) = pd(06,39) - rrt(239) * density(06) 
  pd(34,06) = pd(34,06) + rrt(239) * density(39) 
  pd(34,39) = pd(34,39) + rrt(239) * density(06) 
  pd(39,06) = pd(39,06) - rrt(239) * density(39) 
  pd(39,39) = pd(39,39) - rrt(239) * density(06) 
  pd(03,03) = pd(03,03) - rrt(240) * density(40) 
  pd(03,40) = pd(03,40) - rrt(240) * density(03) 
  pd(35,03) = pd(35,03) + rrt(240) * density(40) 
  pd(35,40) = pd(35,40) + rrt(240) * density(03) 
  pd(40,03) = pd(40,03) - rrt(240) * density(40) 
  pd(40,40) = pd(40,40) - rrt(240) * density(03) 
  pd(04,04) = pd(04,04) - rrt(241) * density(40) 
  pd(04,40) = pd(04,40) - rrt(241) * density(04) 
  pd(35,04) = pd(35,04) + rrt(241) * density(40) 
  pd(35,40) = pd(35,40) + rrt(241) * density(04) 
  pd(40,04) = pd(40,04) - rrt(241) * density(40) 
  pd(40,40) = pd(40,40) - rrt(241) * density(04) 
  pd(05,05) = pd(05,05) - rrt(242) * density(40) 
  pd(05,40) = pd(05,40) - rrt(242) * density(05) 
  pd(35,05) = pd(35,05) + rrt(242) * density(40) 
  pd(35,40) = pd(35,40) + rrt(242) * density(05) 
  pd(40,05) = pd(40,05) - rrt(242) * density(40) 
  pd(40,40) = pd(40,40) - rrt(242) * density(05) 
  pd(06,06) = pd(06,06) - rrt(243) * density(40) 
  pd(06,40) = pd(06,40) - rrt(243) * density(06) 
  pd(35,06) = pd(35,06) + rrt(243) * density(40) 
  pd(35,40) = pd(35,40) + rrt(243) * density(06) 
  pd(40,06) = pd(40,06) - rrt(243) * density(40) 
  pd(40,40) = pd(40,40) - rrt(243) * density(06) 
  pd(03,03) = pd(03,03) - rrt(244) * density(41) 
  pd(03,41) = pd(03,41) - rrt(244) * density(03) 
  pd(15,03) = pd(15,03) + rrt(244) * density(41) 
  pd(15,41) = pd(15,41) + rrt(244) * density(03) 
  pd(34,03) = pd(34,03) + rrt(244) * density(41) 
  pd(34,41) = pd(34,41) + rrt(244) * density(03) 
  pd(41,03) = pd(41,03) - rrt(244) * density(41) 
  pd(41,41) = pd(41,41) - rrt(244) * density(03) 
  pd(04,04) = pd(04,04) - rrt(245) * density(41) 
  pd(04,41) = pd(04,41) - rrt(245) * density(04) 
  pd(15,04) = pd(15,04) + rrt(245) * density(41) 
  pd(15,41) = pd(15,41) + rrt(245) * density(04) 
  pd(34,04) = pd(34,04) + rrt(245) * density(41) 
  pd(34,41) = pd(34,41) + rrt(245) * density(04) 
  pd(41,04) = pd(41,04) - rrt(245) * density(41) 
  pd(41,41) = pd(41,41) - rrt(245) * density(04) 
  pd(05,05) = pd(05,05) - rrt(246) * density(41) 
  pd(05,41) = pd(05,41) - rrt(246) * density(05) 
  pd(15,05) = pd(15,05) + rrt(246) * density(41) 
  pd(15,41) = pd(15,41) + rrt(246) * density(05) 
  pd(34,05) = pd(34,05) + rrt(246) * density(41) 
  pd(34,41) = pd(34,41) + rrt(246) * density(05) 
  pd(41,05) = pd(41,05) - rrt(246) * density(41) 
  pd(41,41) = pd(41,41) - rrt(246) * density(05) 
  pd(06,06) = pd(06,06) - rrt(247) * density(41) 
  pd(06,41) = pd(06,41) - rrt(247) * density(06) 
  pd(15,06) = pd(15,06) + rrt(247) * density(41) 
  pd(15,41) = pd(15,41) + rrt(247) * density(06) 
  pd(34,06) = pd(34,06) + rrt(247) * density(41) 
  pd(34,41) = pd(34,41) + rrt(247) * density(06) 
  pd(41,06) = pd(41,06) - rrt(247) * density(41) 
  pd(41,41) = pd(41,41) - rrt(247) * density(06) 
  pd(03,03) = pd(03,03) - rrt(248) * density(42) 
  pd(03,42) = pd(03,42) - rrt(248) * density(03) 
  pd(15,03) = pd(15,03) + rrt(248) * density(42) 
  pd(15,42) = pd(15,42) + rrt(248) * density(03) 
  pd(35,03) = pd(35,03) + rrt(248) * density(42) 
  pd(35,42) = pd(35,42) + rrt(248) * density(03) 
  pd(42,03) = pd(42,03) - rrt(248) * density(42) 
  pd(42,42) = pd(42,42) - rrt(248) * density(03) 
  pd(04,04) = pd(04,04) - rrt(249) * density(42) 
  pd(04,42) = pd(04,42) - rrt(249) * density(04) 
  pd(15,04) = pd(15,04) + rrt(249) * density(42) 
  pd(15,42) = pd(15,42) + rrt(249) * density(04) 
  pd(35,04) = pd(35,04) + rrt(249) * density(42) 
  pd(35,42) = pd(35,42) + rrt(249) * density(04) 
  pd(42,04) = pd(42,04) - rrt(249) * density(42) 
  pd(42,42) = pd(42,42) - rrt(249) * density(04) 
  pd(05,05) = pd(05,05) - rrt(250) * density(42) 
  pd(05,42) = pd(05,42) - rrt(250) * density(05) 
  pd(15,05) = pd(15,05) + rrt(250) * density(42) 
  pd(15,42) = pd(15,42) + rrt(250) * density(05) 
  pd(35,05) = pd(35,05) + rrt(250) * density(42) 
  pd(35,42) = pd(35,42) + rrt(250) * density(05) 
  pd(42,05) = pd(42,05) - rrt(250) * density(42) 
  pd(42,42) = pd(42,42) - rrt(250) * density(05) 
  pd(06,06) = pd(06,06) - rrt(251) * density(42) 
  pd(06,42) = pd(06,42) - rrt(251) * density(06) 
  pd(15,06) = pd(15,06) + rrt(251) * density(42) 
  pd(15,42) = pd(15,42) + rrt(251) * density(06) 
  pd(35,06) = pd(35,06) + rrt(251) * density(42) 
  pd(35,42) = pd(35,42) + rrt(251) * density(06) 
  pd(42,06) = pd(42,06) - rrt(251) * density(42) 
  pd(42,42) = pd(42,42) - rrt(251) * density(06) 
  pd(03,03) = pd(03,03) - rrt(252) * density(26) 
  pd(03,26) = pd(03,26) - rrt(252) * density(03) 
  pd(26,03) = pd(26,03) - rrt(252) * density(26) 
  pd(26,26) = pd(26,26) - rrt(252) * density(03) 
  pd(31,03) = pd(31,03) + rrt(252) * density(26) 
  pd(31,26) = pd(31,26) + rrt(252) * density(03) 
  pd(04,04) = pd(04,04) - rrt(253) * density(26) 
  pd(04,26) = pd(04,26) - rrt(253) * density(04) 
  pd(26,04) = pd(26,04) - rrt(253) * density(26) 
  pd(26,26) = pd(26,26) - rrt(253) * density(04) 
  pd(31,04) = pd(31,04) + rrt(253) * density(26) 
  pd(31,26) = pd(31,26) + rrt(253) * density(04) 
  pd(05,05) = pd(05,05) - rrt(254) * density(26) 
  pd(05,26) = pd(05,26) - rrt(254) * density(05) 
  pd(26,05) = pd(26,05) - rrt(254) * density(26) 
  pd(26,26) = pd(26,26) - rrt(254) * density(05) 
  pd(31,05) = pd(31,05) + rrt(254) * density(26) 
  pd(31,26) = pd(31,26) + rrt(254) * density(05) 
  pd(06,06) = pd(06,06) - rrt(255) * density(26) 
  pd(06,26) = pd(06,26) - rrt(255) * density(06) 
  pd(26,06) = pd(26,06) - rrt(255) * density(26) 
  pd(26,26) = pd(26,26) - rrt(255) * density(06) 
  pd(31,06) = pd(31,06) + rrt(255) * density(26) 
  pd(31,26) = pd(31,26) + rrt(255) * density(06) 
  pd(03,03) = pd(03,03) - rrt(256) * density(36) 
  pd(03,36) = pd(03,36) - rrt(256) * density(03) 
  pd(27,03) = pd(27,03) + rrt(256) * density(36) 
  pd(27,36) = pd(27,36) + rrt(256) * density(03) 
  pd(36,03) = pd(36,03) - rrt(256) * density(36) 
  pd(36,36) = pd(36,36) - rrt(256) * density(03) 
  pd(04,04) = pd(04,04) - rrt(257) * density(36) 
  pd(04,36) = pd(04,36) - rrt(257) * density(04) 
  pd(27,04) = pd(27,04) + rrt(257) * density(36) 
  pd(27,36) = pd(27,36) + rrt(257) * density(04) 
  pd(36,04) = pd(36,04) - rrt(257) * density(36) 
  pd(36,36) = pd(36,36) - rrt(257) * density(04) 
  pd(05,05) = pd(05,05) - rrt(258) * density(36) 
  pd(05,36) = pd(05,36) - rrt(258) * density(05) 
  pd(27,05) = pd(27,05) + rrt(258) * density(36) 
  pd(27,36) = pd(27,36) + rrt(258) * density(05) 
  pd(36,05) = pd(36,05) - rrt(258) * density(36) 
  pd(36,36) = pd(36,36) - rrt(258) * density(05) 
  pd(06,06) = pd(06,06) - rrt(259) * density(36) 
  pd(06,36) = pd(06,36) - rrt(259) * density(06) 
  pd(27,06) = pd(27,06) + rrt(259) * density(36) 
  pd(27,36) = pd(27,36) + rrt(259) * density(06) 
  pd(36,06) = pd(36,06) - rrt(259) * density(36) 
  pd(36,36) = pd(36,36) - rrt(259) * density(06) 
  pd(07,07) = pd(07,07) - rrt(260) * density(37) 
  pd(07,37) = pd(07,37) - rrt(260) * density(07) 
  pd(09,07) = pd(09,07) + rrt(260) * density(37) 
  pd(09,37) = pd(09,37) + rrt(260) * density(07) 
  pd(36,07) = pd(36,07) + rrt(260) * density(37) 
  pd(36,37) = pd(36,37) + rrt(260) * density(07) 
  pd(37,07) = pd(37,07) - rrt(260) * density(37) 
  pd(37,37) = pd(37,37) - rrt(260) * density(07) 
  pd(08,08) = pd(08,08) - rrt(261) * density(37) 
  pd(08,37) = pd(08,37) - rrt(261) * density(08) 
  pd(09,08) = pd(09,08) + rrt(261) * density(37) 
  pd(09,37) = pd(09,37) + rrt(261) * density(08) 
  pd(36,08) = pd(36,08) + rrt(261) * density(37) 
  pd(36,37) = pd(36,37) + rrt(261) * density(08) 
  pd(37,08) = pd(37,08) - rrt(261) * density(37) 
  pd(37,37) = pd(37,37) - rrt(261) * density(08) 
  pd(07,07) = pd(07,07) - rrt(262) * density(37) 
  pd(07,37) = pd(07,37) - rrt(262) * density(07) 
  pd(25,07) = pd(25,07) + rrt(262) * density(37) 
  pd(25,37) = pd(25,37) + rrt(262) * density(07) 
  pd(26,07) = pd(26,07) + rrt(262) * density(37) 
  pd(26,37) = pd(26,37) + rrt(262) * density(07) 
  pd(37,07) = pd(37,07) - rrt(262) * density(37) 
  pd(37,37) = pd(37,37) - rrt(262) * density(07) 
  pd(08,08) = pd(08,08) - rrt(263) * density(37) 
  pd(08,37) = pd(08,37) - rrt(263) * density(08) 
  pd(25,08) = pd(25,08) + rrt(263) * density(37) 
  pd(25,37) = pd(25,37) + rrt(263) * density(08) 
  pd(26,08) = pd(26,08) + rrt(263) * density(37) 
  pd(26,37) = pd(26,37) + rrt(263) * density(08) 
  pd(37,08) = pd(37,08) - rrt(263) * density(37) 
  pd(37,37) = pd(37,37) - rrt(263) * density(08) 
  pd(07,07) = pd(07,07) - rrt(264) * density(32) 
  pd(07,32) = pd(07,32) - rrt(264) * density(07) 
  pd(09,07) = pd(09,07) + rrt(264) * density(32) 
  pd(09,32) = pd(09,32) + rrt(264) * density(07) 
  pd(28,07) = pd(28,07) + rrt(264) * density(32) 
  pd(28,32) = pd(28,32) + rrt(264) * density(07) 
  pd(32,07) = pd(32,07) - rrt(264) * density(32) 
  pd(32,32) = pd(32,32) - rrt(264) * density(07) 
  pd(08,08) = pd(08,08) - rrt(265) * density(32) 
  pd(08,32) = pd(08,32) - rrt(265) * density(08) 
  pd(09,08) = pd(09,08) + rrt(265) * density(32) 
  pd(09,32) = pd(09,32) + rrt(265) * density(08) 
  pd(28,08) = pd(28,08) + rrt(265) * density(32) 
  pd(28,32) = pd(28,32) + rrt(265) * density(08) 
  pd(32,08) = pd(32,08) - rrt(265) * density(32) 
  pd(32,32) = pd(32,32) - rrt(265) * density(08) 
  pd(07,07) = pd(07,07) - rrt(266) * density(28) 
  pd(07,28) = pd(07,28) - rrt(266) * density(07) 
  pd(09,07) = pd(09,07) + rrt(266) * density(28) 
  pd(09,28) = pd(09,28) + rrt(266) * density(07) 
  pd(26,07) = pd(26,07) + rrt(266) * density(28) 
  pd(26,28) = pd(26,28) + rrt(266) * density(07) 
  pd(28,07) = pd(28,07) - rrt(266) * density(28) 
  pd(28,28) = pd(28,28) - rrt(266) * density(07) 
  pd(08,08) = pd(08,08) - rrt(267) * density(28) 
  pd(08,28) = pd(08,28) - rrt(267) * density(08) 
  pd(09,08) = pd(09,08) + rrt(267) * density(28) 
  pd(09,28) = pd(09,28) + rrt(267) * density(08) 
  pd(26,08) = pd(26,08) + rrt(267) * density(28) 
  pd(26,28) = pd(26,28) + rrt(267) * density(08) 
  pd(28,08) = pd(28,08) - rrt(267) * density(28) 
  pd(28,28) = pd(28,28) - rrt(267) * density(08) 
  pd(07,07) = pd(07,07) - rrt(268) * density(39) 
  pd(07,39) = pd(07,39) - rrt(268) * density(07) 
  pd(34,07) = pd(34,07) + rrt(268) * density(39) 
  pd(34,39) = pd(34,39) + rrt(268) * density(07) 
  pd(39,07) = pd(39,07) - rrt(268) * density(39) 
  pd(39,39) = pd(39,39) - rrt(268) * density(07) 
  pd(08,08) = pd(08,08) - rrt(269) * density(39) 
  pd(08,39) = pd(08,39) - rrt(269) * density(08) 
  pd(34,08) = pd(34,08) + rrt(269) * density(39) 
  pd(34,39) = pd(34,39) + rrt(269) * density(08) 
  pd(39,08) = pd(39,08) - rrt(269) * density(39) 
  pd(39,39) = pd(39,39) - rrt(269) * density(08) 
  pd(07,07) = pd(07,07) - rrt(270) * density(40) 
  pd(07,40) = pd(07,40) - rrt(270) * density(07) 
  pd(35,07) = pd(35,07) + rrt(270) * density(40) 
  pd(35,40) = pd(35,40) + rrt(270) * density(07) 
  pd(40,07) = pd(40,07) - rrt(270) * density(40) 
  pd(40,40) = pd(40,40) - rrt(270) * density(07) 
  pd(08,08) = pd(08,08) - rrt(271) * density(40) 
  pd(08,40) = pd(08,40) - rrt(271) * density(08) 
  pd(35,08) = pd(35,08) + rrt(271) * density(40) 
  pd(35,40) = pd(35,40) + rrt(271) * density(08) 
  pd(40,08) = pd(40,08) - rrt(271) * density(40) 
  pd(40,40) = pd(40,40) - rrt(271) * density(08) 
  pd(07,07) = pd(07,07) - rrt(272) * density(41) 
  pd(07,41) = pd(07,41) - rrt(272) * density(07) 
  pd(15,07) = pd(15,07) + rrt(272) * density(41) 
  pd(15,41) = pd(15,41) + rrt(272) * density(07) 
  pd(34,07) = pd(34,07) + rrt(272) * density(41) 
  pd(34,41) = pd(34,41) + rrt(272) * density(07) 
  pd(41,07) = pd(41,07) - rrt(272) * density(41) 
  pd(41,41) = pd(41,41) - rrt(272) * density(07) 
  pd(08,08) = pd(08,08) - rrt(273) * density(41) 
  pd(08,41) = pd(08,41) - rrt(273) * density(08) 
  pd(15,08) = pd(15,08) + rrt(273) * density(41) 
  pd(15,41) = pd(15,41) + rrt(273) * density(08) 
  pd(34,08) = pd(34,08) + rrt(273) * density(41) 
  pd(34,41) = pd(34,41) + rrt(273) * density(08) 
  pd(41,08) = pd(41,08) - rrt(273) * density(41) 
  pd(41,41) = pd(41,41) - rrt(273) * density(08) 
  pd(07,07) = pd(07,07) - rrt(274) * density(42) 
  pd(07,42) = pd(07,42) - rrt(274) * density(07) 
  pd(15,07) = pd(15,07) + rrt(274) * density(42) 
  pd(15,42) = pd(15,42) + rrt(274) * density(07) 
  pd(35,07) = pd(35,07) + rrt(274) * density(42) 
  pd(35,42) = pd(35,42) + rrt(274) * density(07) 
  pd(42,07) = pd(42,07) - rrt(274) * density(42) 
  pd(42,42) = pd(42,42) - rrt(274) * density(07) 
  pd(08,08) = pd(08,08) - rrt(275) * density(42) 
  pd(08,42) = pd(08,42) - rrt(275) * density(08) 
  pd(15,08) = pd(15,08) + rrt(275) * density(42) 
  pd(15,42) = pd(15,42) + rrt(275) * density(08) 
  pd(35,08) = pd(35,08) + rrt(275) * density(42) 
  pd(35,42) = pd(35,42) + rrt(275) * density(08) 
  pd(42,08) = pd(42,08) - rrt(275) * density(42) 
  pd(42,42) = pd(42,42) - rrt(275) * density(08) 
  pd(07,07) = pd(07,07) - rrt(276) * density(26) 
  pd(07,26) = pd(07,26) - rrt(276) * density(07) 
  pd(26,07) = pd(26,07) - rrt(276) * density(26) 
  pd(26,26) = pd(26,26) - rrt(276) * density(07) 
  pd(31,07) = pd(31,07) + rrt(276) * density(26) 
  pd(31,26) = pd(31,26) + rrt(276) * density(07) 
  pd(08,08) = pd(08,08) - rrt(277) * density(26) 
  pd(08,26) = pd(08,26) - rrt(277) * density(08) 
  pd(26,08) = pd(26,08) - rrt(277) * density(26) 
  pd(26,26) = pd(26,26) - rrt(277) * density(08) 
  pd(31,08) = pd(31,08) + rrt(277) * density(26) 
  pd(31,26) = pd(31,26) + rrt(277) * density(08) 
  pd(07,07) = pd(07,07) - rrt(278) * density(36) 
  pd(07,36) = pd(07,36) - rrt(278) * density(07) 
  pd(27,07) = pd(27,07) + rrt(278) * density(36) 
  pd(27,36) = pd(27,36) + rrt(278) * density(07) 
  pd(36,07) = pd(36,07) - rrt(278) * density(36) 
  pd(36,36) = pd(36,36) - rrt(278) * density(07) 
  pd(08,08) = pd(08,08) - rrt(279) * density(36) 
  pd(08,36) = pd(08,36) - rrt(279) * density(08) 
  pd(27,08) = pd(27,08) + rrt(279) * density(36) 
  pd(27,36) = pd(27,36) + rrt(279) * density(08) 
  pd(36,08) = pd(36,08) - rrt(279) * density(36) 
  pd(36,36) = pd(36,36) - rrt(279) * density(08) 
  pd(07,07) = pd(07,07) - rrt(280) * density(36) 
  pd(07,36) = pd(07,36) - rrt(280) * density(07) 
  pd(27,07) = pd(27,07) + rrt(280) * density(36) 
  pd(27,36) = pd(27,36) + rrt(280) * density(07) 
  pd(36,07) = pd(36,07) - rrt(280) * density(36) 
  pd(36,36) = pd(36,36) - rrt(280) * density(07) 
  pd(08,08) = pd(08,08) - rrt(281) * density(36) 
  pd(08,36) = pd(08,36) - rrt(281) * density(08) 
  pd(27,08) = pd(27,08) + rrt(281) * density(36) 
  pd(27,36) = pd(27,36) + rrt(281) * density(08) 
  pd(36,08) = pd(36,08) - rrt(281) * density(36) 
  pd(36,36) = pd(36,36) - rrt(281) * density(08) 
  pd(09,09) = pd(09,09) - rrt(282) * density(37) 
  pd(09,37) = pd(09,37) - rrt(282) * density(09) 
  pd(25,09) = pd(25,09) + rrt(282) * density(37) 
  pd(25,37) = pd(25,37) + rrt(282) * density(09) 
  pd(28,09) = pd(28,09) + rrt(282) * density(37) 
  pd(28,37) = pd(28,37) + rrt(282) * density(09) 
  pd(37,09) = pd(37,09) - rrt(282) * density(37) 
  pd(37,37) = pd(37,37) - rrt(282) * density(09) 
  pd(01,09) = pd(01,09) + rrt(283) * density(39) 
  pd(01,39) = pd(01,39) + rrt(283) * density(09) 
  pd(02,09) = pd(02,09) + rrt(283) * density(39) 
  pd(02,39) = pd(02,39) + rrt(283) * density(09) 
  pd(09,09) = pd(09,09) - rrt(283) * density(39) 
  pd(09,39) = pd(09,39) - rrt(283) * density(09) 
  pd(39,09) = pd(39,09) - rrt(283) * density(39) 
  pd(39,39) = pd(39,39) - rrt(283) * density(09) 
  pd(01,09) = pd(01,09) + rrt(284) * density(34) 
  pd(01,34) = pd(01,34) + rrt(284) * density(09) 
  pd(02,09) = pd(02,09) + rrt(284) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(284) * density(09) * 2.0d0
  pd(09,09) = pd(09,09) - rrt(284) * density(34) 
  pd(09,34) = pd(09,34) - rrt(284) * density(09) 
  pd(34,09) = pd(34,09) - rrt(284) * density(34) 
  pd(34,34) = pd(34,34) - rrt(284) * density(09) 
  pd(02,09) = pd(02,09) + rrt(285) * density(30) 
  pd(02,30) = pd(02,30) + rrt(285) * density(09) 
  pd(09,09) = pd(09,09) - rrt(285) * density(30) 
  pd(09,30) = pd(09,30) - rrt(285) * density(09) 
  pd(29,09) = pd(29,09) + rrt(285) * density(30) 
  pd(29,30) = pd(29,30) + rrt(285) * density(09) 
  pd(30,09) = pd(30,09) - rrt(285) * density(30) 
  pd(30,30) = pd(30,30) - rrt(285) * density(09) 
  pd(02,09) = pd(02,09) + rrt(286) * density(31) 
  pd(02,31) = pd(02,31) + rrt(286) * density(09) 
  pd(09,09) = pd(09,09) - rrt(286) * density(31) 
  pd(09,31) = pd(09,31) - rrt(286) * density(09) 
  pd(30,09) = pd(30,09) + rrt(286) * density(31) 
  pd(30,31) = pd(30,31) + rrt(286) * density(09) 
  pd(31,09) = pd(31,09) - rrt(286) * density(31) 
  pd(31,31) = pd(31,31) - rrt(286) * density(09) 
  pd(02,09) = pd(02,09) + rrt(287) * density(30) 
  pd(02,30) = pd(02,30) + rrt(287) * density(09) 
  pd(09,09) = pd(09,09) - rrt(287) * density(30) 
  pd(09,30) = pd(09,30) - rrt(287) * density(09) 
  pd(29,09) = pd(29,09) + rrt(287) * density(30) 
  pd(29,30) = pd(29,30) + rrt(287) * density(09) 
  pd(30,09) = pd(30,09) - rrt(287) * density(30) 
  pd(30,30) = pd(30,30) - rrt(287) * density(09) 
  pd(02,09) = pd(02,09) + rrt(288) * density(31) 
  pd(02,31) = pd(02,31) + rrt(288) * density(09) 
  pd(09,09) = pd(09,09) - rrt(288) * density(31) 
  pd(09,31) = pd(09,31) - rrt(288) * density(09) 
  pd(30,09) = pd(30,09) + rrt(288) * density(31) 
  pd(30,31) = pd(30,31) + rrt(288) * density(09) 
  pd(31,09) = pd(31,09) - rrt(288) * density(31) 
  pd(31,31) = pd(31,31) - rrt(288) * density(09) 
  pd(09,09) = pd(09,09) - rrt(289) * density(32) 
  pd(09,32) = pd(09,32) - rrt(289) * density(09) 
  pd(22,09) = pd(22,09) + rrt(289) * density(32) 
  pd(22,32) = pd(22,32) + rrt(289) * density(09) 
  pd(28,09) = pd(28,09) + rrt(289) * density(32) 
  pd(28,32) = pd(28,32) + rrt(289) * density(09) 
  pd(32,09) = pd(32,09) - rrt(289) * density(32) 
  pd(32,32) = pd(32,32) - rrt(289) * density(09) 
  pd(10,10) = pd(10,10) - rrt(290) * density(37) 
  pd(10,37) = pd(10,37) - rrt(290) * density(10) 
  pd(25,10) = pd(25,10) + rrt(290) * density(37) 
  pd(25,37) = pd(25,37) + rrt(290) * density(10) 
  pd(28,10) = pd(28,10) + rrt(290) * density(37) 
  pd(28,37) = pd(28,37) + rrt(290) * density(10) 
  pd(37,10) = pd(37,10) - rrt(290) * density(37) 
  pd(37,37) = pd(37,37) - rrt(290) * density(10) 
  pd(01,10) = pd(01,10) + rrt(291) * density(39) 
  pd(01,39) = pd(01,39) + rrt(291) * density(10) 
  pd(02,10) = pd(02,10) + rrt(291) * density(39) 
  pd(02,39) = pd(02,39) + rrt(291) * density(10) 
  pd(10,10) = pd(10,10) - rrt(291) * density(39) 
  pd(10,39) = pd(10,39) - rrt(291) * density(10) 
  pd(39,10) = pd(39,10) - rrt(291) * density(39) 
  pd(39,39) = pd(39,39) - rrt(291) * density(10) 
  pd(01,10) = pd(01,10) + rrt(292) * density(34) 
  pd(01,34) = pd(01,34) + rrt(292) * density(10) 
  pd(02,10) = pd(02,10) + rrt(292) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(292) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(292) * density(34) 
  pd(10,34) = pd(10,34) - rrt(292) * density(10) 
  pd(34,10) = pd(34,10) - rrt(292) * density(34) 
  pd(34,34) = pd(34,34) - rrt(292) * density(10) 
  pd(02,10) = pd(02,10) + rrt(293) * density(30) 
  pd(02,30) = pd(02,30) + rrt(293) * density(10) 
  pd(10,10) = pd(10,10) - rrt(293) * density(30) 
  pd(10,30) = pd(10,30) - rrt(293) * density(10) 
  pd(29,10) = pd(29,10) + rrt(293) * density(30) 
  pd(29,30) = pd(29,30) + rrt(293) * density(10) 
  pd(30,10) = pd(30,10) - rrt(293) * density(30) 
  pd(30,30) = pd(30,30) - rrt(293) * density(10) 
  pd(02,10) = pd(02,10) + rrt(294) * density(31) 
  pd(02,31) = pd(02,31) + rrt(294) * density(10) 
  pd(10,10) = pd(10,10) - rrt(294) * density(31) 
  pd(10,31) = pd(10,31) - rrt(294) * density(10) 
  pd(30,10) = pd(30,10) + rrt(294) * density(31) 
  pd(30,31) = pd(30,31) + rrt(294) * density(10) 
  pd(31,10) = pd(31,10) - rrt(294) * density(31) 
  pd(31,31) = pd(31,31) - rrt(294) * density(10) 
  pd(02,10) = pd(02,10) + rrt(295) * density(30) 
  pd(02,30) = pd(02,30) + rrt(295) * density(10) 
  pd(10,10) = pd(10,10) - rrt(295) * density(30) 
  pd(10,30) = pd(10,30) - rrt(295) * density(10) 
  pd(29,10) = pd(29,10) + rrt(295) * density(30) 
  pd(29,30) = pd(29,30) + rrt(295) * density(10) 
  pd(30,10) = pd(30,10) - rrt(295) * density(30) 
  pd(30,30) = pd(30,30) - rrt(295) * density(10) 
  pd(02,10) = pd(02,10) + rrt(296) * density(31) 
  pd(02,31) = pd(02,31) + rrt(296) * density(10) 
  pd(10,10) = pd(10,10) - rrt(296) * density(31) 
  pd(10,31) = pd(10,31) - rrt(296) * density(10) 
  pd(30,10) = pd(30,10) + rrt(296) * density(31) 
  pd(30,31) = pd(30,31) + rrt(296) * density(10) 
  pd(31,10) = pd(31,10) - rrt(296) * density(31) 
  pd(31,31) = pd(31,31) - rrt(296) * density(10) 
  pd(10,10) = pd(10,10) - rrt(297) * density(32) 
  pd(10,32) = pd(10,32) - rrt(297) * density(10) 
  pd(22,10) = pd(22,10) + rrt(297) * density(32) 
  pd(22,32) = pd(22,32) + rrt(297) * density(10) 
  pd(28,10) = pd(28,10) + rrt(297) * density(32) 
  pd(28,32) = pd(28,32) + rrt(297) * density(10) 
  pd(32,10) = pd(32,10) - rrt(297) * density(32) 
  pd(32,32) = pd(32,32) - rrt(297) * density(10) 
  pd(01,11) = pd(01,11) + rrt(298) * density(39) 
  pd(01,39) = pd(01,39) + rrt(298) * density(11) 
  pd(02,11) = pd(02,11) + rrt(298) * density(39) 
  pd(02,39) = pd(02,39) + rrt(298) * density(11) 
  pd(11,11) = pd(11,11) - rrt(298) * density(39) 
  pd(11,39) = pd(11,39) - rrt(298) * density(11) 
  pd(39,11) = pd(39,11) - rrt(298) * density(39) 
  pd(39,39) = pd(39,39) - rrt(298) * density(11) 
  pd(01,12) = pd(01,12) + rrt(299) * density(39) 
  pd(01,39) = pd(01,39) + rrt(299) * density(12) 
  pd(02,12) = pd(02,12) + rrt(299) * density(39) 
  pd(02,39) = pd(02,39) + rrt(299) * density(12) 
  pd(12,12) = pd(12,12) - rrt(299) * density(39) 
  pd(12,39) = pd(12,39) - rrt(299) * density(12) 
  pd(39,12) = pd(39,12) - rrt(299) * density(39) 
  pd(39,39) = pd(39,39) - rrt(299) * density(12) 
  pd(01,13) = pd(01,13) + rrt(300) * density(39) 
  pd(01,39) = pd(01,39) + rrt(300) * density(13) 
  pd(02,13) = pd(02,13) + rrt(300) * density(39) 
  pd(02,39) = pd(02,39) + rrt(300) * density(13) 
  pd(13,13) = pd(13,13) - rrt(300) * density(39) 
  pd(13,39) = pd(13,39) - rrt(300) * density(13) 
  pd(39,13) = pd(39,13) - rrt(300) * density(39) 
  pd(39,39) = pd(39,39) - rrt(300) * density(13) 
  pd(01,14) = pd(01,14) + rrt(301) * density(39) 
  pd(01,39) = pd(01,39) + rrt(301) * density(14) 
  pd(02,14) = pd(02,14) + rrt(301) * density(39) 
  pd(02,39) = pd(02,39) + rrt(301) * density(14) 
  pd(14,14) = pd(14,14) - rrt(301) * density(39) 
  pd(14,39) = pd(14,39) - rrt(301) * density(14) 
  pd(39,14) = pd(39,14) - rrt(301) * density(39) 
  pd(39,39) = pd(39,39) - rrt(301) * density(14) 
  pd(01,11) = pd(01,11) + rrt(302) * density(34) 
  pd(01,34) = pd(01,34) + rrt(302) * density(11) 
  pd(02,11) = pd(02,11) + rrt(302) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(302) * density(11) * 2.0d0
  pd(11,11) = pd(11,11) - rrt(302) * density(34) 
  pd(11,34) = pd(11,34) - rrt(302) * density(11) 
  pd(34,11) = pd(34,11) - rrt(302) * density(34) 
  pd(34,34) = pd(34,34) - rrt(302) * density(11) 
  pd(01,12) = pd(01,12) + rrt(303) * density(34) 
  pd(01,34) = pd(01,34) + rrt(303) * density(12) 
  pd(02,12) = pd(02,12) + rrt(303) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(303) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(303) * density(34) 
  pd(12,34) = pd(12,34) - rrt(303) * density(12) 
  pd(34,12) = pd(34,12) - rrt(303) * density(34) 
  pd(34,34) = pd(34,34) - rrt(303) * density(12) 
  pd(01,13) = pd(01,13) + rrt(304) * density(34) 
  pd(01,34) = pd(01,34) + rrt(304) * density(13) 
  pd(02,13) = pd(02,13) + rrt(304) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(304) * density(13) * 2.0d0
  pd(13,13) = pd(13,13) - rrt(304) * density(34) 
  pd(13,34) = pd(13,34) - rrt(304) * density(13) 
  pd(34,13) = pd(34,13) - rrt(304) * density(34) 
  pd(34,34) = pd(34,34) - rrt(304) * density(13) 
  pd(01,14) = pd(01,14) + rrt(305) * density(34) 
  pd(01,34) = pd(01,34) + rrt(305) * density(14) 
  pd(02,14) = pd(02,14) + rrt(305) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(305) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(305) * density(34) 
  pd(14,34) = pd(14,34) - rrt(305) * density(14) 
  pd(34,14) = pd(34,14) - rrt(305) * density(34) 
  pd(34,34) = pd(34,34) - rrt(305) * density(14) 
  pd(02,11) = pd(02,11) + rrt(306) * density(30) 
  pd(02,30) = pd(02,30) + rrt(306) * density(11) 
  pd(11,11) = pd(11,11) - rrt(306) * density(30) 
  pd(11,30) = pd(11,30) - rrt(306) * density(11) 
  pd(29,11) = pd(29,11) + rrt(306) * density(30) 
  pd(29,30) = pd(29,30) + rrt(306) * density(11) 
  pd(30,11) = pd(30,11) - rrt(306) * density(30) 
  pd(30,30) = pd(30,30) - rrt(306) * density(11) 
  pd(02,12) = pd(02,12) + rrt(307) * density(30) 
  pd(02,30) = pd(02,30) + rrt(307) * density(12) 
  pd(12,12) = pd(12,12) - rrt(307) * density(30) 
  pd(12,30) = pd(12,30) - rrt(307) * density(12) 
  pd(29,12) = pd(29,12) + rrt(307) * density(30) 
  pd(29,30) = pd(29,30) + rrt(307) * density(12) 
  pd(30,12) = pd(30,12) - rrt(307) * density(30) 
  pd(30,30) = pd(30,30) - rrt(307) * density(12) 
  pd(02,13) = pd(02,13) + rrt(308) * density(30) 
  pd(02,30) = pd(02,30) + rrt(308) * density(13) 
  pd(13,13) = pd(13,13) - rrt(308) * density(30) 
  pd(13,30) = pd(13,30) - rrt(308) * density(13) 
  pd(29,13) = pd(29,13) + rrt(308) * density(30) 
  pd(29,30) = pd(29,30) + rrt(308) * density(13) 
  pd(30,13) = pd(30,13) - rrt(308) * density(30) 
  pd(30,30) = pd(30,30) - rrt(308) * density(13) 
  pd(02,14) = pd(02,14) + rrt(309) * density(30) 
  pd(02,30) = pd(02,30) + rrt(309) * density(14) 
  pd(14,14) = pd(14,14) - rrt(309) * density(30) 
  pd(14,30) = pd(14,30) - rrt(309) * density(14) 
  pd(29,14) = pd(29,14) + rrt(309) * density(30) 
  pd(29,30) = pd(29,30) + rrt(309) * density(14) 
  pd(30,14) = pd(30,14) - rrt(309) * density(30) 
  pd(30,30) = pd(30,30) - rrt(309) * density(14) 
  pd(02,11) = pd(02,11) + rrt(310) * density(31) 
  pd(02,31) = pd(02,31) + rrt(310) * density(11) 
  pd(11,11) = pd(11,11) - rrt(310) * density(31) 
  pd(11,31) = pd(11,31) - rrt(310) * density(11) 
  pd(30,11) = pd(30,11) + rrt(310) * density(31) 
  pd(30,31) = pd(30,31) + rrt(310) * density(11) 
  pd(31,11) = pd(31,11) - rrt(310) * density(31) 
  pd(31,31) = pd(31,31) - rrt(310) * density(11) 
  pd(02,12) = pd(02,12) + rrt(311) * density(31) 
  pd(02,31) = pd(02,31) + rrt(311) * density(12) 
  pd(12,12) = pd(12,12) - rrt(311) * density(31) 
  pd(12,31) = pd(12,31) - rrt(311) * density(12) 
  pd(30,12) = pd(30,12) + rrt(311) * density(31) 
  pd(30,31) = pd(30,31) + rrt(311) * density(12) 
  pd(31,12) = pd(31,12) - rrt(311) * density(31) 
  pd(31,31) = pd(31,31) - rrt(311) * density(12) 
  pd(02,13) = pd(02,13) + rrt(312) * density(31) 
  pd(02,31) = pd(02,31) + rrt(312) * density(13) 
  pd(13,13) = pd(13,13) - rrt(312) * density(31) 
  pd(13,31) = pd(13,31) - rrt(312) * density(13) 
  pd(30,13) = pd(30,13) + rrt(312) * density(31) 
  pd(30,31) = pd(30,31) + rrt(312) * density(13) 
  pd(31,13) = pd(31,13) - rrt(312) * density(31) 
  pd(31,31) = pd(31,31) - rrt(312) * density(13) 
  pd(02,14) = pd(02,14) + rrt(313) * density(31) 
  pd(02,31) = pd(02,31) + rrt(313) * density(14) 
  pd(14,14) = pd(14,14) - rrt(313) * density(31) 
  pd(14,31) = pd(14,31) - rrt(313) * density(14) 
  pd(30,14) = pd(30,14) + rrt(313) * density(31) 
  pd(30,31) = pd(30,31) + rrt(313) * density(14) 
  pd(31,14) = pd(31,14) - rrt(313) * density(31) 
  pd(31,31) = pd(31,31) - rrt(313) * density(14) 
  pd(02,11) = pd(02,11) + rrt(314) * density(30) 
  pd(02,30) = pd(02,30) + rrt(314) * density(11) 
  pd(11,11) = pd(11,11) - rrt(314) * density(30) 
  pd(11,30) = pd(11,30) - rrt(314) * density(11) 
  pd(29,11) = pd(29,11) + rrt(314) * density(30) 
  pd(29,30) = pd(29,30) + rrt(314) * density(11) 
  pd(30,11) = pd(30,11) - rrt(314) * density(30) 
  pd(30,30) = pd(30,30) - rrt(314) * density(11) 
  pd(02,12) = pd(02,12) + rrt(315) * density(30) 
  pd(02,30) = pd(02,30) + rrt(315) * density(12) 
  pd(12,12) = pd(12,12) - rrt(315) * density(30) 
  pd(12,30) = pd(12,30) - rrt(315) * density(12) 
  pd(29,12) = pd(29,12) + rrt(315) * density(30) 
  pd(29,30) = pd(29,30) + rrt(315) * density(12) 
  pd(30,12) = pd(30,12) - rrt(315) * density(30) 
  pd(30,30) = pd(30,30) - rrt(315) * density(12) 
  pd(02,13) = pd(02,13) + rrt(316) * density(30) 
  pd(02,30) = pd(02,30) + rrt(316) * density(13) 
  pd(13,13) = pd(13,13) - rrt(316) * density(30) 
  pd(13,30) = pd(13,30) - rrt(316) * density(13) 
  pd(29,13) = pd(29,13) + rrt(316) * density(30) 
  pd(29,30) = pd(29,30) + rrt(316) * density(13) 
  pd(30,13) = pd(30,13) - rrt(316) * density(30) 
  pd(30,30) = pd(30,30) - rrt(316) * density(13) 
  pd(02,14) = pd(02,14) + rrt(317) * density(30) 
  pd(02,30) = pd(02,30) + rrt(317) * density(14) 
  pd(14,14) = pd(14,14) - rrt(317) * density(30) 
  pd(14,30) = pd(14,30) - rrt(317) * density(14) 
  pd(29,14) = pd(29,14) + rrt(317) * density(30) 
  pd(29,30) = pd(29,30) + rrt(317) * density(14) 
  pd(30,14) = pd(30,14) - rrt(317) * density(30) 
  pd(30,30) = pd(30,30) - rrt(317) * density(14) 
  pd(02,11) = pd(02,11) + rrt(318) * density(31) 
  pd(02,31) = pd(02,31) + rrt(318) * density(11) 
  pd(11,11) = pd(11,11) - rrt(318) * density(31) 
  pd(11,31) = pd(11,31) - rrt(318) * density(11) 
  pd(30,11) = pd(30,11) + rrt(318) * density(31) 
  pd(30,31) = pd(30,31) + rrt(318) * density(11) 
  pd(31,11) = pd(31,11) - rrt(318) * density(31) 
  pd(31,31) = pd(31,31) - rrt(318) * density(11) 
  pd(02,12) = pd(02,12) + rrt(319) * density(31) 
  pd(02,31) = pd(02,31) + rrt(319) * density(12) 
  pd(12,12) = pd(12,12) - rrt(319) * density(31) 
  pd(12,31) = pd(12,31) - rrt(319) * density(12) 
  pd(30,12) = pd(30,12) + rrt(319) * density(31) 
  pd(30,31) = pd(30,31) + rrt(319) * density(12) 
  pd(31,12) = pd(31,12) - rrt(319) * density(31) 
  pd(31,31) = pd(31,31) - rrt(319) * density(12) 
  pd(02,13) = pd(02,13) + rrt(320) * density(31) 
  pd(02,31) = pd(02,31) + rrt(320) * density(13) 
  pd(13,13) = pd(13,13) - rrt(320) * density(31) 
  pd(13,31) = pd(13,31) - rrt(320) * density(13) 
  pd(30,13) = pd(30,13) + rrt(320) * density(31) 
  pd(30,31) = pd(30,31) + rrt(320) * density(13) 
  pd(31,13) = pd(31,13) - rrt(320) * density(31) 
  pd(31,31) = pd(31,31) - rrt(320) * density(13) 
  pd(02,14) = pd(02,14) + rrt(321) * density(31) 
  pd(02,31) = pd(02,31) + rrt(321) * density(14) 
  pd(14,14) = pd(14,14) - rrt(321) * density(31) 
  pd(14,31) = pd(14,31) - rrt(321) * density(14) 
  pd(30,14) = pd(30,14) + rrt(321) * density(31) 
  pd(30,31) = pd(30,31) + rrt(321) * density(14) 
  pd(31,14) = pd(31,14) - rrt(321) * density(31) 
  pd(31,31) = pd(31,31) - rrt(321) * density(14) 
  pd(11,11) = pd(11,11) - rrt(322) * density(32) 
  pd(11,32) = pd(11,32) - rrt(322) * density(11) 
  pd(22,11) = pd(22,11) + rrt(322) * density(32) 
  pd(22,32) = pd(22,32) + rrt(322) * density(11) 
  pd(28,11) = pd(28,11) + rrt(322) * density(32) 
  pd(28,32) = pd(28,32) + rrt(322) * density(11) 
  pd(32,11) = pd(32,11) - rrt(322) * density(32) 
  pd(32,32) = pd(32,32) - rrt(322) * density(11) 
  pd(12,12) = pd(12,12) - rrt(323) * density(32) 
  pd(12,32) = pd(12,32) - rrt(323) * density(12) 
  pd(22,12) = pd(22,12) + rrt(323) * density(32) 
  pd(22,32) = pd(22,32) + rrt(323) * density(12) 
  pd(28,12) = pd(28,12) + rrt(323) * density(32) 
  pd(28,32) = pd(28,32) + rrt(323) * density(12) 
  pd(32,12) = pd(32,12) - rrt(323) * density(32) 
  pd(32,32) = pd(32,32) - rrt(323) * density(12) 
  pd(13,13) = pd(13,13) - rrt(324) * density(32) 
  pd(13,32) = pd(13,32) - rrt(324) * density(13) 
  pd(22,13) = pd(22,13) + rrt(324) * density(32) 
  pd(22,32) = pd(22,32) + rrt(324) * density(13) 
  pd(28,13) = pd(28,13) + rrt(324) * density(32) 
  pd(28,32) = pd(28,32) + rrt(324) * density(13) 
  pd(32,13) = pd(32,13) - rrt(324) * density(32) 
  pd(32,32) = pd(32,32) - rrt(324) * density(13) 
  pd(14,14) = pd(14,14) - rrt(325) * density(32) 
  pd(14,32) = pd(14,32) - rrt(325) * density(14) 
  pd(22,14) = pd(22,14) + rrt(325) * density(32) 
  pd(22,32) = pd(22,32) + rrt(325) * density(14) 
  pd(28,14) = pd(28,14) + rrt(325) * density(32) 
  pd(28,32) = pd(28,32) + rrt(325) * density(14) 
  pd(32,14) = pd(32,14) - rrt(325) * density(32) 
  pd(32,32) = pd(32,32) - rrt(325) * density(14) 
  pd(02,15) = pd(02,15) + rrt(326) * density(26) 
  pd(02,26) = pd(02,26) + rrt(326) * density(15) 
  pd(15,15) = pd(15,15) - rrt(326) * density(26) 
  pd(15,26) = pd(15,26) - rrt(326) * density(15) 
  pd(26,15) = pd(26,15) - rrt(326) * density(26) 
  pd(26,26) = pd(26,26) - rrt(326) * density(15) 
  pd(36,15) = pd(36,15) + rrt(326) * density(26) 
  pd(36,26) = pd(36,26) + rrt(326) * density(15) 
  pd(09,15) = pd(09,15) + rrt(327) * density(28) 
  pd(09,28) = pd(09,28) + rrt(327) * density(15) 
  pd(15,15) = pd(15,15) - rrt(327) * density(28) 
  pd(15,28) = pd(15,28) - rrt(327) * density(15) 
  pd(28,15) = pd(28,15) - rrt(327) * density(28) 
  pd(28,28) = pd(28,28) - rrt(327) * density(15) 
  pd(36,15) = pd(36,15) + rrt(327) * density(28) 
  pd(36,28) = pd(36,28) + rrt(327) * density(15) 
  pd(09,15) = pd(09,15) + rrt(328) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(328) * density(15) * 2.0d0
  pd(15,15) = pd(15,15) - rrt(328) * density(29) 
  pd(15,29) = pd(15,29) - rrt(328) * density(15) 
  pd(29,15) = pd(29,15) - rrt(328) * density(29) 
  pd(29,29) = pd(29,29) - rrt(328) * density(15) 
  pd(36,15) = pd(36,15) + rrt(328) * density(29) 
  pd(36,29) = pd(36,29) + rrt(328) * density(15) 
  pd(09,15) = pd(09,15) + rrt(329) * density(32) 
  pd(09,32) = pd(09,32) + rrt(329) * density(15) 
  pd(15,15) = pd(15,15) - rrt(329) * density(32) 
  pd(15,32) = pd(15,32) - rrt(329) * density(15) 
  pd(32,15) = pd(32,15) - rrt(329) * density(32) 
  pd(32,32) = pd(32,32) - rrt(329) * density(15) 
  pd(37,15) = pd(37,15) + rrt(329) * density(32) 
  pd(37,32) = pd(37,32) + rrt(329) * density(15) 
  pd(15,15) = pd(15,15) - rrt(330) * density(32) 
  pd(15,32) = pd(15,32) - rrt(330) * density(15) 
  pd(25,15) = pd(25,15) + rrt(330) * density(32) 
  pd(25,32) = pd(25,32) + rrt(330) * density(15) 
  pd(28,15) = pd(28,15) + rrt(330) * density(32) 
  pd(28,32) = pd(28,32) + rrt(330) * density(15) 
  pd(32,15) = pd(32,15) - rrt(330) * density(32) 
  pd(32,32) = pd(32,32) - rrt(330) * density(15) 
  pd(15,15) = pd(15,15) - rrt(331) * density(37) 
  pd(15,37) = pd(15,37) - rrt(331) * density(15) 
  pd(25,15) = pd(25,15) + rrt(331) * density(37) 
  pd(25,37) = pd(25,37) + rrt(331) * density(15) 
  pd(36,15) = pd(36,15) + rrt(331) * density(37) 
  pd(36,37) = pd(36,37) + rrt(331) * density(15) 
  pd(37,15) = pd(37,15) - rrt(331) * density(37) 
  pd(37,37) = pd(37,37) - rrt(331) * density(15) 
  pd(15,15) = pd(15,15) - rrt(332) * density(36) 
  pd(15,36) = pd(15,36) - rrt(332) * density(15) 
  pd(36,15) = pd(36,15) - rrt(332) * density(36) 
  pd(36,36) = pd(36,36) - rrt(332) * density(15) 
  pd(38,15) = pd(38,15) + rrt(332) * density(36) 
  pd(38,36) = pd(38,36) + rrt(332) * density(15) 
  pd(15,15) = pd(15,15) - rrt(333) * density(40) 
  pd(15,40) = pd(15,40) - rrt(333) * density(15) 
  pd(40,15) = pd(40,15) - rrt(333) * density(40) 
  pd(40,40) = pd(40,40) - rrt(333) * density(15) 
  pd(42,15) = pd(42,15) + rrt(333) * density(40) 
  pd(42,40) = pd(42,40) + rrt(333) * density(15) 
  pd(01,15) = pd(01,15) + rrt(334) * density(39) 
  pd(01,39) = pd(01,39) + rrt(334) * density(15) 
  pd(15,15) = pd(15,15) - rrt(334) * density(39) 
  pd(15,39) = pd(15,39) - rrt(334) * density(15) 
  pd(24,15) = pd(24,15) + rrt(334) * density(39) 
  pd(24,39) = pd(24,39) + rrt(334) * density(15) 
  pd(39,15) = pd(39,15) - rrt(334) * density(39) 
  pd(39,39) = pd(39,39) - rrt(334) * density(15) 
  pd(15,15) = pd(15,15) - rrt(335) * density(39) 
  pd(15,39) = pd(15,39) - rrt(335) * density(15) 
  pd(39,15) = pd(39,15) - rrt(335) * density(39) 
  pd(39,39) = pd(39,39) - rrt(335) * density(15) 
  pd(41,15) = pd(41,15) + rrt(335) * density(39) 
  pd(41,39) = pd(41,39) + rrt(335) * density(15) 
  pd(01,15) = pd(01,15) + rrt(336) * density(40) 
  pd(01,40) = pd(01,40) + rrt(336) * density(15) 
  pd(15,15) = pd(15,15) + rrt(336) * density(40) 
  pd(15,40) = pd(15,40) + rrt(336) * density(15) 
  pd(40,15) = pd(40,15) - rrt(336) * density(40) 
  pd(40,40) = pd(40,40) - rrt(336) * density(15) 
  pd(01,15) = pd(01,15) + rrt(337) * density(41) 
  pd(01,41) = pd(01,41) + rrt(337) * density(15) 
  pd(24,15) = pd(24,15) + rrt(337) * density(41) 
  pd(24,41) = pd(24,41) + rrt(337) * density(15) 
  pd(41,15) = pd(41,15) - rrt(337) * density(41) 
  pd(41,41) = pd(41,41) - rrt(337) * density(15) 
  pd(02,16) = pd(02,16) + rrt(338) * density(26) 
  pd(02,26) = pd(02,26) + rrt(338) * density(16) 
  pd(16,16) = pd(16,16) - rrt(338) * density(26) 
  pd(16,26) = pd(16,26) - rrt(338) * density(16) 
  pd(26,16) = pd(26,16) - rrt(338) * density(26) 
  pd(26,26) = pd(26,26) - rrt(338) * density(16) 
  pd(36,16) = pd(36,16) + rrt(338) * density(26) 
  pd(36,26) = pd(36,26) + rrt(338) * density(16) 
  pd(02,17) = pd(02,17) + rrt(339) * density(26) 
  pd(02,26) = pd(02,26) + rrt(339) * density(17) 
  pd(17,17) = pd(17,17) - rrt(339) * density(26) 
  pd(17,26) = pd(17,26) - rrt(339) * density(17) 
  pd(26,17) = pd(26,17) - rrt(339) * density(26) 
  pd(26,26) = pd(26,26) - rrt(339) * density(17) 
  pd(36,17) = pd(36,17) + rrt(339) * density(26) 
  pd(36,26) = pd(36,26) + rrt(339) * density(17) 
  pd(02,18) = pd(02,18) + rrt(340) * density(26) 
  pd(02,26) = pd(02,26) + rrt(340) * density(18) 
  pd(18,18) = pd(18,18) - rrt(340) * density(26) 
  pd(18,26) = pd(18,26) - rrt(340) * density(18) 
  pd(26,18) = pd(26,18) - rrt(340) * density(26) 
  pd(26,26) = pd(26,26) - rrt(340) * density(18) 
  pd(36,18) = pd(36,18) + rrt(340) * density(26) 
  pd(36,26) = pd(36,26) + rrt(340) * density(18) 
  pd(09,16) = pd(09,16) + rrt(341) * density(28) 
  pd(09,28) = pd(09,28) + rrt(341) * density(16) 
  pd(16,16) = pd(16,16) - rrt(341) * density(28) 
  pd(16,28) = pd(16,28) - rrt(341) * density(16) 
  pd(28,16) = pd(28,16) - rrt(341) * density(28) 
  pd(28,28) = pd(28,28) - rrt(341) * density(16) 
  pd(36,16) = pd(36,16) + rrt(341) * density(28) 
  pd(36,28) = pd(36,28) + rrt(341) * density(16) 
  pd(09,17) = pd(09,17) + rrt(342) * density(28) 
  pd(09,28) = pd(09,28) + rrt(342) * density(17) 
  pd(17,17) = pd(17,17) - rrt(342) * density(28) 
  pd(17,28) = pd(17,28) - rrt(342) * density(17) 
  pd(28,17) = pd(28,17) - rrt(342) * density(28) 
  pd(28,28) = pd(28,28) - rrt(342) * density(17) 
  pd(36,17) = pd(36,17) + rrt(342) * density(28) 
  pd(36,28) = pd(36,28) + rrt(342) * density(17) 
  pd(09,18) = pd(09,18) + rrt(343) * density(28) 
  pd(09,28) = pd(09,28) + rrt(343) * density(18) 
  pd(18,18) = pd(18,18) - rrt(343) * density(28) 
  pd(18,28) = pd(18,28) - rrt(343) * density(18) 
  pd(28,18) = pd(28,18) - rrt(343) * density(28) 
  pd(28,28) = pd(28,28) - rrt(343) * density(18) 
  pd(36,18) = pd(36,18) + rrt(343) * density(28) 
  pd(36,28) = pd(36,28) + rrt(343) * density(18) 
  pd(09,16) = pd(09,16) + rrt(344) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(344) * density(16) * 2.0d0
  pd(16,16) = pd(16,16) - rrt(344) * density(29) 
  pd(16,29) = pd(16,29) - rrt(344) * density(16) 
  pd(29,16) = pd(29,16) - rrt(344) * density(29) 
  pd(29,29) = pd(29,29) - rrt(344) * density(16) 
  pd(36,16) = pd(36,16) + rrt(344) * density(29) 
  pd(36,29) = pd(36,29) + rrt(344) * density(16) 
  pd(09,17) = pd(09,17) + rrt(345) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(345) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(345) * density(29) 
  pd(17,29) = pd(17,29) - rrt(345) * density(17) 
  pd(29,17) = pd(29,17) - rrt(345) * density(29) 
  pd(29,29) = pd(29,29) - rrt(345) * density(17) 
  pd(36,17) = pd(36,17) + rrt(345) * density(29) 
  pd(36,29) = pd(36,29) + rrt(345) * density(17) 
  pd(09,18) = pd(09,18) + rrt(346) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(346) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(346) * density(29) 
  pd(18,29) = pd(18,29) - rrt(346) * density(18) 
  pd(29,18) = pd(29,18) - rrt(346) * density(29) 
  pd(29,29) = pd(29,29) - rrt(346) * density(18) 
  pd(36,18) = pd(36,18) + rrt(346) * density(29) 
  pd(36,29) = pd(36,29) + rrt(346) * density(18) 
  pd(09,16) = pd(09,16) + rrt(347) * density(32) 
  pd(09,32) = pd(09,32) + rrt(347) * density(16) 
  pd(16,16) = pd(16,16) - rrt(347) * density(32) 
  pd(16,32) = pd(16,32) - rrt(347) * density(16) 
  pd(32,16) = pd(32,16) - rrt(347) * density(32) 
  pd(32,32) = pd(32,32) - rrt(347) * density(16) 
  pd(37,16) = pd(37,16) + rrt(347) * density(32) 
  pd(37,32) = pd(37,32) + rrt(347) * density(16) 
  pd(09,17) = pd(09,17) + rrt(348) * density(32) 
  pd(09,32) = pd(09,32) + rrt(348) * density(17) 
  pd(17,17) = pd(17,17) - rrt(348) * density(32) 
  pd(17,32) = pd(17,32) - rrt(348) * density(17) 
  pd(32,17) = pd(32,17) - rrt(348) * density(32) 
  pd(32,32) = pd(32,32) - rrt(348) * density(17) 
  pd(37,17) = pd(37,17) + rrt(348) * density(32) 
  pd(37,32) = pd(37,32) + rrt(348) * density(17) 
  pd(09,18) = pd(09,18) + rrt(349) * density(32) 
  pd(09,32) = pd(09,32) + rrt(349) * density(18) 
  pd(18,18) = pd(18,18) - rrt(349) * density(32) 
  pd(18,32) = pd(18,32) - rrt(349) * density(18) 
  pd(32,18) = pd(32,18) - rrt(349) * density(32) 
  pd(32,32) = pd(32,32) - rrt(349) * density(18) 
  pd(37,18) = pd(37,18) + rrt(349) * density(32) 
  pd(37,32) = pd(37,32) + rrt(349) * density(18) 
  pd(16,16) = pd(16,16) - rrt(350) * density(32) 
  pd(16,32) = pd(16,32) - rrt(350) * density(16) 
  pd(25,16) = pd(25,16) + rrt(350) * density(32) 
  pd(25,32) = pd(25,32) + rrt(350) * density(16) 
  pd(28,16) = pd(28,16) + rrt(350) * density(32) 
  pd(28,32) = pd(28,32) + rrt(350) * density(16) 
  pd(32,16) = pd(32,16) - rrt(350) * density(32) 
  pd(32,32) = pd(32,32) - rrt(350) * density(16) 
  pd(17,17) = pd(17,17) - rrt(351) * density(32) 
  pd(17,32) = pd(17,32) - rrt(351) * density(17) 
  pd(25,17) = pd(25,17) + rrt(351) * density(32) 
  pd(25,32) = pd(25,32) + rrt(351) * density(17) 
  pd(28,17) = pd(28,17) + rrt(351) * density(32) 
  pd(28,32) = pd(28,32) + rrt(351) * density(17) 
  pd(32,17) = pd(32,17) - rrt(351) * density(32) 
  pd(32,32) = pd(32,32) - rrt(351) * density(17) 
  pd(18,18) = pd(18,18) - rrt(352) * density(32) 
  pd(18,32) = pd(18,32) - rrt(352) * density(18) 
  pd(25,18) = pd(25,18) + rrt(352) * density(32) 
  pd(25,32) = pd(25,32) + rrt(352) * density(18) 
  pd(28,18) = pd(28,18) + rrt(352) * density(32) 
  pd(28,32) = pd(28,32) + rrt(352) * density(18) 
  pd(32,18) = pd(32,18) - rrt(352) * density(32) 
  pd(32,32) = pd(32,32) - rrt(352) * density(18) 
  pd(16,16) = pd(16,16) - rrt(353) * density(37) 
  pd(16,37) = pd(16,37) - rrt(353) * density(16) 
  pd(25,16) = pd(25,16) + rrt(353) * density(37) 
  pd(25,37) = pd(25,37) + rrt(353) * density(16) 
  pd(36,16) = pd(36,16) + rrt(353) * density(37) 
  pd(36,37) = pd(36,37) + rrt(353) * density(16) 
  pd(37,16) = pd(37,16) - rrt(353) * density(37) 
  pd(37,37) = pd(37,37) - rrt(353) * density(16) 
  pd(17,17) = pd(17,17) - rrt(354) * density(37) 
  pd(17,37) = pd(17,37) - rrt(354) * density(17) 
  pd(25,17) = pd(25,17) + rrt(354) * density(37) 
  pd(25,37) = pd(25,37) + rrt(354) * density(17) 
  pd(36,17) = pd(36,17) + rrt(354) * density(37) 
  pd(36,37) = pd(36,37) + rrt(354) * density(17) 
  pd(37,17) = pd(37,17) - rrt(354) * density(37) 
  pd(37,37) = pd(37,37) - rrt(354) * density(17) 
  pd(18,18) = pd(18,18) - rrt(355) * density(37) 
  pd(18,37) = pd(18,37) - rrt(355) * density(18) 
  pd(25,18) = pd(25,18) + rrt(355) * density(37) 
  pd(25,37) = pd(25,37) + rrt(355) * density(18) 
  pd(36,18) = pd(36,18) + rrt(355) * density(37) 
  pd(36,37) = pd(36,37) + rrt(355) * density(18) 
  pd(37,18) = pd(37,18) - rrt(355) * density(37) 
  pd(37,37) = pd(37,37) - rrt(355) * density(18) 
  pd(16,16) = pd(16,16) - rrt(356) * density(36) 
  pd(16,36) = pd(16,36) - rrt(356) * density(16) 
  pd(36,16) = pd(36,16) - rrt(356) * density(36) 
  pd(36,36) = pd(36,36) - rrt(356) * density(16) 
  pd(38,16) = pd(38,16) + rrt(356) * density(36) 
  pd(38,36) = pd(38,36) + rrt(356) * density(16) 
  pd(17,17) = pd(17,17) - rrt(357) * density(36) 
  pd(17,36) = pd(17,36) - rrt(357) * density(17) 
  pd(36,17) = pd(36,17) - rrt(357) * density(36) 
  pd(36,36) = pd(36,36) - rrt(357) * density(17) 
  pd(38,17) = pd(38,17) + rrt(357) * density(36) 
  pd(38,36) = pd(38,36) + rrt(357) * density(17) 
  pd(18,18) = pd(18,18) - rrt(358) * density(36) 
  pd(18,36) = pd(18,36) - rrt(358) * density(18) 
  pd(36,18) = pd(36,18) - rrt(358) * density(36) 
  pd(36,36) = pd(36,36) - rrt(358) * density(18) 
  pd(38,18) = pd(38,18) + rrt(358) * density(36) 
  pd(38,36) = pd(38,36) + rrt(358) * density(18) 
  pd(16,16) = pd(16,16) - rrt(359) * density(40) 
  pd(16,40) = pd(16,40) - rrt(359) * density(16) 
  pd(40,16) = pd(40,16) - rrt(359) * density(40) 
  pd(40,40) = pd(40,40) - rrt(359) * density(16) 
  pd(42,16) = pd(42,16) + rrt(359) * density(40) 
  pd(42,40) = pd(42,40) + rrt(359) * density(16) 
  pd(17,17) = pd(17,17) - rrt(360) * density(40) 
  pd(17,40) = pd(17,40) - rrt(360) * density(17) 
  pd(40,17) = pd(40,17) - rrt(360) * density(40) 
  pd(40,40) = pd(40,40) - rrt(360) * density(17) 
  pd(42,17) = pd(42,17) + rrt(360) * density(40) 
  pd(42,40) = pd(42,40) + rrt(360) * density(17) 
  pd(18,18) = pd(18,18) - rrt(361) * density(40) 
  pd(18,40) = pd(18,40) - rrt(361) * density(18) 
  pd(40,18) = pd(40,18) - rrt(361) * density(40) 
  pd(40,40) = pd(40,40) - rrt(361) * density(18) 
  pd(42,18) = pd(42,18) + rrt(361) * density(40) 
  pd(42,40) = pd(42,40) + rrt(361) * density(18) 
  pd(01,16) = pd(01,16) + rrt(362) * density(39) 
  pd(01,39) = pd(01,39) + rrt(362) * density(16) 
  pd(16,16) = pd(16,16) - rrt(362) * density(39) 
  pd(16,39) = pd(16,39) - rrt(362) * density(16) 
  pd(24,16) = pd(24,16) + rrt(362) * density(39) 
  pd(24,39) = pd(24,39) + rrt(362) * density(16) 
  pd(39,16) = pd(39,16) - rrt(362) * density(39) 
  pd(39,39) = pd(39,39) - rrt(362) * density(16) 
  pd(01,17) = pd(01,17) + rrt(363) * density(39) 
  pd(01,39) = pd(01,39) + rrt(363) * density(17) 
  pd(17,17) = pd(17,17) - rrt(363) * density(39) 
  pd(17,39) = pd(17,39) - rrt(363) * density(17) 
  pd(24,17) = pd(24,17) + rrt(363) * density(39) 
  pd(24,39) = pd(24,39) + rrt(363) * density(17) 
  pd(39,17) = pd(39,17) - rrt(363) * density(39) 
  pd(39,39) = pd(39,39) - rrt(363) * density(17) 
  pd(01,18) = pd(01,18) + rrt(364) * density(39) 
  pd(01,39) = pd(01,39) + rrt(364) * density(18) 
  pd(18,18) = pd(18,18) - rrt(364) * density(39) 
  pd(18,39) = pd(18,39) - rrt(364) * density(18) 
  pd(24,18) = pd(24,18) + rrt(364) * density(39) 
  pd(24,39) = pd(24,39) + rrt(364) * density(18) 
  pd(39,18) = pd(39,18) - rrt(364) * density(39) 
  pd(39,39) = pd(39,39) - rrt(364) * density(18) 
  pd(16,16) = pd(16,16) - rrt(365) * density(39) 
  pd(16,39) = pd(16,39) - rrt(365) * density(16) 
  pd(39,16) = pd(39,16) - rrt(365) * density(39) 
  pd(39,39) = pd(39,39) - rrt(365) * density(16) 
  pd(41,16) = pd(41,16) + rrt(365) * density(39) 
  pd(41,39) = pd(41,39) + rrt(365) * density(16) 
  pd(17,17) = pd(17,17) - rrt(366) * density(39) 
  pd(17,39) = pd(17,39) - rrt(366) * density(17) 
  pd(39,17) = pd(39,17) - rrt(366) * density(39) 
  pd(39,39) = pd(39,39) - rrt(366) * density(17) 
  pd(41,17) = pd(41,17) + rrt(366) * density(39) 
  pd(41,39) = pd(41,39) + rrt(366) * density(17) 
  pd(18,18) = pd(18,18) - rrt(367) * density(39) 
  pd(18,39) = pd(18,39) - rrt(367) * density(18) 
  pd(39,18) = pd(39,18) - rrt(367) * density(39) 
  pd(39,39) = pd(39,39) - rrt(367) * density(18) 
  pd(41,18) = pd(41,18) + rrt(367) * density(39) 
  pd(41,39) = pd(41,39) + rrt(367) * density(18) 
  pd(01,16) = pd(01,16) + rrt(368) * density(40) 
  pd(01,40) = pd(01,40) + rrt(368) * density(16) 
  pd(15,16) = pd(15,16) + rrt(368) * density(40) 
  pd(15,40) = pd(15,40) + rrt(368) * density(16) 
  pd(40,16) = pd(40,16) - rrt(368) * density(40) 
  pd(40,40) = pd(40,40) - rrt(368) * density(16) 
  pd(01,17) = pd(01,17) + rrt(369) * density(40) 
  pd(01,40) = pd(01,40) + rrt(369) * density(17) 
  pd(15,17) = pd(15,17) + rrt(369) * density(40) 
  pd(15,40) = pd(15,40) + rrt(369) * density(17) 
  pd(40,17) = pd(40,17) - rrt(369) * density(40) 
  pd(40,40) = pd(40,40) - rrt(369) * density(17) 
  pd(01,18) = pd(01,18) + rrt(370) * density(40) 
  pd(01,40) = pd(01,40) + rrt(370) * density(18) 
  pd(15,18) = pd(15,18) + rrt(370) * density(40) 
  pd(15,40) = pd(15,40) + rrt(370) * density(18) 
  pd(40,18) = pd(40,18) - rrt(370) * density(40) 
  pd(40,40) = pd(40,40) - rrt(370) * density(18) 
  pd(01,16) = pd(01,16) + rrt(371) * density(41) 
  pd(01,41) = pd(01,41) + rrt(371) * density(16) 
  pd(24,16) = pd(24,16) + rrt(371) * density(41) 
  pd(24,41) = pd(24,41) + rrt(371) * density(16) 
  pd(41,16) = pd(41,16) - rrt(371) * density(41) 
  pd(41,41) = pd(41,41) - rrt(371) * density(16) 
  pd(01,17) = pd(01,17) + rrt(372) * density(41) 
  pd(01,41) = pd(01,41) + rrt(372) * density(17) 
  pd(24,17) = pd(24,17) + rrt(372) * density(41) 
  pd(24,41) = pd(24,41) + rrt(372) * density(17) 
  pd(41,17) = pd(41,17) - rrt(372) * density(41) 
  pd(41,41) = pd(41,41) - rrt(372) * density(17) 
  pd(01,18) = pd(01,18) + rrt(373) * density(41) 
  pd(01,41) = pd(01,41) + rrt(373) * density(18) 
  pd(24,18) = pd(24,18) + rrt(373) * density(41) 
  pd(24,41) = pd(24,41) + rrt(373) * density(18) 
  pd(41,18) = pd(41,18) - rrt(373) * density(41) 
  pd(41,41) = pd(41,41) - rrt(373) * density(18) 
  pd(02,19) = pd(02,19) + rrt(374) * density(26) 
  pd(02,26) = pd(02,26) + rrt(374) * density(19) 
  pd(19,19) = pd(19,19) - rrt(374) * density(26) 
  pd(19,26) = pd(19,26) - rrt(374) * density(19) 
  pd(26,19) = pd(26,19) - rrt(374) * density(26) 
  pd(26,26) = pd(26,26) - rrt(374) * density(19) 
  pd(36,19) = pd(36,19) + rrt(374) * density(26) 
  pd(36,26) = pd(36,26) + rrt(374) * density(19) 
  pd(02,20) = pd(02,20) + rrt(375) * density(26) 
  pd(02,26) = pd(02,26) + rrt(375) * density(20) 
  pd(20,20) = pd(20,20) - rrt(375) * density(26) 
  pd(20,26) = pd(20,26) - rrt(375) * density(20) 
  pd(26,20) = pd(26,20) - rrt(375) * density(26) 
  pd(26,26) = pd(26,26) - rrt(375) * density(20) 
  pd(36,20) = pd(36,20) + rrt(375) * density(26) 
  pd(36,26) = pd(36,26) + rrt(375) * density(20) 
  pd(09,19) = pd(09,19) + rrt(376) * density(28) 
  pd(09,28) = pd(09,28) + rrt(376) * density(19) 
  pd(19,19) = pd(19,19) - rrt(376) * density(28) 
  pd(19,28) = pd(19,28) - rrt(376) * density(19) 
  pd(28,19) = pd(28,19) - rrt(376) * density(28) 
  pd(28,28) = pd(28,28) - rrt(376) * density(19) 
  pd(36,19) = pd(36,19) + rrt(376) * density(28) 
  pd(36,28) = pd(36,28) + rrt(376) * density(19) 
  pd(09,20) = pd(09,20) + rrt(377) * density(28) 
  pd(09,28) = pd(09,28) + rrt(377) * density(20) 
  pd(20,20) = pd(20,20) - rrt(377) * density(28) 
  pd(20,28) = pd(20,28) - rrt(377) * density(20) 
  pd(28,20) = pd(28,20) - rrt(377) * density(28) 
  pd(28,28) = pd(28,28) - rrt(377) * density(20) 
  pd(36,20) = pd(36,20) + rrt(377) * density(28) 
  pd(36,28) = pd(36,28) + rrt(377) * density(20) 
  pd(09,19) = pd(09,19) + rrt(378) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(378) * density(19) * 2.0d0
  pd(19,19) = pd(19,19) - rrt(378) * density(29) 
  pd(19,29) = pd(19,29) - rrt(378) * density(19) 
  pd(29,19) = pd(29,19) - rrt(378) * density(29) 
  pd(29,29) = pd(29,29) - rrt(378) * density(19) 
  pd(36,19) = pd(36,19) + rrt(378) * density(29) 
  pd(36,29) = pd(36,29) + rrt(378) * density(19) 
  pd(09,20) = pd(09,20) + rrt(379) * density(29) * 2.0d0
  pd(09,29) = pd(09,29) + rrt(379) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(379) * density(29) 
  pd(20,29) = pd(20,29) - rrt(379) * density(20) 
  pd(29,20) = pd(29,20) - rrt(379) * density(29) 
  pd(29,29) = pd(29,29) - rrt(379) * density(20) 
  pd(36,20) = pd(36,20) + rrt(379) * density(29) 
  pd(36,29) = pd(36,29) + rrt(379) * density(20) 
  pd(09,19) = pd(09,19) + rrt(380) * density(32) 
  pd(09,32) = pd(09,32) + rrt(380) * density(19) 
  pd(19,19) = pd(19,19) - rrt(380) * density(32) 
  pd(19,32) = pd(19,32) - rrt(380) * density(19) 
  pd(32,19) = pd(32,19) - rrt(380) * density(32) 
  pd(32,32) = pd(32,32) - rrt(380) * density(19) 
  pd(37,19) = pd(37,19) + rrt(380) * density(32) 
  pd(37,32) = pd(37,32) + rrt(380) * density(19) 
  pd(09,20) = pd(09,20) + rrt(381) * density(32) 
  pd(09,32) = pd(09,32) + rrt(381) * density(20) 
  pd(20,20) = pd(20,20) - rrt(381) * density(32) 
  pd(20,32) = pd(20,32) - rrt(381) * density(20) 
  pd(32,20) = pd(32,20) - rrt(381) * density(32) 
  pd(32,32) = pd(32,32) - rrt(381) * density(20) 
  pd(37,20) = pd(37,20) + rrt(381) * density(32) 
  pd(37,32) = pd(37,32) + rrt(381) * density(20) 
  pd(19,19) = pd(19,19) - rrt(382) * density(32) 
  pd(19,32) = pd(19,32) - rrt(382) * density(19) 
  pd(25,19) = pd(25,19) + rrt(382) * density(32) 
  pd(25,32) = pd(25,32) + rrt(382) * density(19) 
  pd(28,19) = pd(28,19) + rrt(382) * density(32) 
  pd(28,32) = pd(28,32) + rrt(382) * density(19) 
  pd(32,19) = pd(32,19) - rrt(382) * density(32) 
  pd(32,32) = pd(32,32) - rrt(382) * density(19) 
  pd(20,20) = pd(20,20) - rrt(383) * density(32) 
  pd(20,32) = pd(20,32) - rrt(383) * density(20) 
  pd(25,20) = pd(25,20) + rrt(383) * density(32) 
  pd(25,32) = pd(25,32) + rrt(383) * density(20) 
  pd(28,20) = pd(28,20) + rrt(383) * density(32) 
  pd(28,32) = pd(28,32) + rrt(383) * density(20) 
  pd(32,20) = pd(32,20) - rrt(383) * density(32) 
  pd(32,32) = pd(32,32) - rrt(383) * density(20) 
  pd(19,19) = pd(19,19) - rrt(384) * density(37) 
  pd(19,37) = pd(19,37) - rrt(384) * density(19) 
  pd(25,19) = pd(25,19) + rrt(384) * density(37) 
  pd(25,37) = pd(25,37) + rrt(384) * density(19) 
  pd(36,19) = pd(36,19) + rrt(384) * density(37) 
  pd(36,37) = pd(36,37) + rrt(384) * density(19) 
  pd(37,19) = pd(37,19) - rrt(384) * density(37) 
  pd(37,37) = pd(37,37) - rrt(384) * density(19) 
  pd(20,20) = pd(20,20) - rrt(385) * density(37) 
  pd(20,37) = pd(20,37) - rrt(385) * density(20) 
  pd(25,20) = pd(25,20) + rrt(385) * density(37) 
  pd(25,37) = pd(25,37) + rrt(385) * density(20) 
  pd(36,20) = pd(36,20) + rrt(385) * density(37) 
  pd(36,37) = pd(36,37) + rrt(385) * density(20) 
  pd(37,20) = pd(37,20) - rrt(385) * density(37) 
  pd(37,37) = pd(37,37) - rrt(385) * density(20) 
  pd(19,19) = pd(19,19) - rrt(386) * density(36) 
  pd(19,36) = pd(19,36) - rrt(386) * density(19) 
  pd(36,19) = pd(36,19) - rrt(386) * density(36) 
  pd(36,36) = pd(36,36) - rrt(386) * density(19) 
  pd(38,19) = pd(38,19) + rrt(386) * density(36) 
  pd(38,36) = pd(38,36) + rrt(386) * density(19) 
  pd(20,20) = pd(20,20) - rrt(387) * density(36) 
  pd(20,36) = pd(20,36) - rrt(387) * density(20) 
  pd(36,20) = pd(36,20) - rrt(387) * density(36) 
  pd(36,36) = pd(36,36) - rrt(387) * density(20) 
  pd(38,20) = pd(38,20) + rrt(387) * density(36) 
  pd(38,36) = pd(38,36) + rrt(387) * density(20) 
  pd(19,19) = pd(19,19) - rrt(388) * density(40) 
  pd(19,40) = pd(19,40) - rrt(388) * density(19) 
  pd(40,19) = pd(40,19) - rrt(388) * density(40) 
  pd(40,40) = pd(40,40) - rrt(388) * density(19) 
  pd(42,19) = pd(42,19) + rrt(388) * density(40) 
  pd(42,40) = pd(42,40) + rrt(388) * density(19) 
  pd(20,20) = pd(20,20) - rrt(389) * density(40) 
  pd(20,40) = pd(20,40) - rrt(389) * density(20) 
  pd(40,20) = pd(40,20) - rrt(389) * density(40) 
  pd(40,40) = pd(40,40) - rrt(389) * density(20) 
  pd(42,20) = pd(42,20) + rrt(389) * density(40) 
  pd(42,40) = pd(42,40) + rrt(389) * density(20) 
  pd(01,19) = pd(01,19) + rrt(390) * density(39) 
  pd(01,39) = pd(01,39) + rrt(390) * density(19) 
  pd(19,19) = pd(19,19) - rrt(390) * density(39) 
  pd(19,39) = pd(19,39) - rrt(390) * density(19) 
  pd(24,19) = pd(24,19) + rrt(390) * density(39) 
  pd(24,39) = pd(24,39) + rrt(390) * density(19) 
  pd(39,19) = pd(39,19) - rrt(390) * density(39) 
  pd(39,39) = pd(39,39) - rrt(390) * density(19) 
  pd(01,20) = pd(01,20) + rrt(391) * density(39) 
  pd(01,39) = pd(01,39) + rrt(391) * density(20) 
  pd(20,20) = pd(20,20) - rrt(391) * density(39) 
  pd(20,39) = pd(20,39) - rrt(391) * density(20) 
  pd(24,20) = pd(24,20) + rrt(391) * density(39) 
  pd(24,39) = pd(24,39) + rrt(391) * density(20) 
  pd(39,20) = pd(39,20) - rrt(391) * density(39) 
  pd(39,39) = pd(39,39) - rrt(391) * density(20) 
  pd(19,19) = pd(19,19) - rrt(392) * density(39) 
  pd(19,39) = pd(19,39) - rrt(392) * density(19) 
  pd(39,19) = pd(39,19) - rrt(392) * density(39) 
  pd(39,39) = pd(39,39) - rrt(392) * density(19) 
  pd(41,19) = pd(41,19) + rrt(392) * density(39) 
  pd(41,39) = pd(41,39) + rrt(392) * density(19) 
  pd(20,20) = pd(20,20) - rrt(393) * density(39) 
  pd(20,39) = pd(20,39) - rrt(393) * density(20) 
  pd(39,20) = pd(39,20) - rrt(393) * density(39) 
  pd(39,39) = pd(39,39) - rrt(393) * density(20) 
  pd(41,20) = pd(41,20) + rrt(393) * density(39) 
  pd(41,39) = pd(41,39) + rrt(393) * density(20) 
  pd(01,19) = pd(01,19) + rrt(394) * density(40) 
  pd(01,40) = pd(01,40) + rrt(394) * density(19) 
  pd(15,19) = pd(15,19) + rrt(394) * density(40) 
  pd(15,40) = pd(15,40) + rrt(394) * density(19) 
  pd(40,19) = pd(40,19) - rrt(394) * density(40) 
  pd(40,40) = pd(40,40) - rrt(394) * density(19) 
  pd(01,20) = pd(01,20) + rrt(395) * density(40) 
  pd(01,40) = pd(01,40) + rrt(395) * density(20) 
  pd(15,20) = pd(15,20) + rrt(395) * density(40) 
  pd(15,40) = pd(15,40) + rrt(395) * density(20) 
  pd(40,20) = pd(40,20) - rrt(395) * density(40) 
  pd(40,40) = pd(40,40) - rrt(395) * density(20) 
  pd(01,19) = pd(01,19) + rrt(396) * density(41) 
  pd(01,41) = pd(01,41) + rrt(396) * density(19) 
  pd(24,19) = pd(24,19) + rrt(396) * density(41) 
  pd(24,41) = pd(24,41) + rrt(396) * density(19) 
  pd(41,19) = pd(41,19) - rrt(396) * density(41) 
  pd(41,41) = pd(41,41) - rrt(396) * density(19) 
  pd(01,20) = pd(01,20) + rrt(397) * density(41) 
  pd(01,41) = pd(01,41) + rrt(397) * density(20) 
  pd(24,20) = pd(24,20) + rrt(397) * density(41) 
  pd(24,41) = pd(24,41) + rrt(397) * density(20) 
  pd(41,20) = pd(41,20) - rrt(397) * density(41) 
  pd(41,41) = pd(41,41) - rrt(397) * density(20) 
  pd(09,22) = pd(09,22) + rrt(398) * density(28) 
  pd(09,28) = pd(09,28) + rrt(398) * density(22) 
  pd(22,22) = pd(22,22) - rrt(398) * density(28) 
  pd(22,28) = pd(22,28) - rrt(398) * density(22) 
  pd(28,22) = pd(28,22) - rrt(398) * density(28) 
  pd(28,28) = pd(28,28) - rrt(398) * density(22) 
  pd(32,22) = pd(32,22) + rrt(398) * density(28) 
  pd(32,28) = pd(32,28) + rrt(398) * density(22) 
  pd(22,22) = pd(22,22) - rrt(399) * density(36) 
  pd(22,36) = pd(22,36) - rrt(399) * density(22) 
  pd(25,22) = pd(25,22) + rrt(399) * density(36) 
  pd(25,36) = pd(25,36) + rrt(399) * density(22) 
  pd(28,22) = pd(28,22) + rrt(399) * density(36) 
  pd(28,36) = pd(28,36) + rrt(399) * density(22) 
  pd(36,22) = pd(36,22) - rrt(399) * density(36) 
  pd(36,36) = pd(36,36) - rrt(399) * density(22) 
  pd(15,22) = pd(15,22) + rrt(400) * density(36) 
  pd(15,36) = pd(15,36) + rrt(400) * density(22) 
  pd(22,22) = pd(22,22) - rrt(400) * density(36) 
  pd(22,36) = pd(22,36) - rrt(400) * density(22) 
  pd(32,22) = pd(32,22) + rrt(400) * density(36) 
  pd(32,36) = pd(32,36) + rrt(400) * density(22) 
  pd(36,22) = pd(36,22) - rrt(400) * density(36) 
  pd(36,36) = pd(36,36) - rrt(400) * density(22) 
  pd(22,22) = pd(22,22) - rrt(401) * density(33) 
  pd(22,33) = pd(22,33) - rrt(401) * density(22) 
  pd(23,22) = pd(23,22) + rrt(401) * density(33) 
  pd(23,33) = pd(23,33) + rrt(401) * density(22) 
  pd(32,22) = pd(32,22) + rrt(401) * density(33) 
  pd(32,33) = pd(32,33) + rrt(401) * density(22) 
  pd(33,22) = pd(33,22) - rrt(401) * density(33) 
  pd(33,33) = pd(33,33) - rrt(401) * density(22) 
  pd(09,25) = pd(09,25) + rrt(402) * density(26) 
  pd(09,26) = pd(09,26) + rrt(402) * density(25) 
  pd(25,25) = pd(25,25) - rrt(402) * density(26) 
  pd(25,26) = pd(25,26) - rrt(402) * density(25) 
  pd(26,25) = pd(26,25) - rrt(402) * density(26) 
  pd(26,26) = pd(26,26) - rrt(402) * density(25) 
  pd(36,25) = pd(36,25) + rrt(402) * density(26) 
  pd(36,26) = pd(36,26) + rrt(402) * density(25) 
  pd(02,25) = pd(02,25) + rrt(403) * density(26) 
  pd(02,26) = pd(02,26) + rrt(403) * density(25) 
  pd(25,25) = pd(25,25) - rrt(403) * density(26) 
  pd(25,26) = pd(25,26) - rrt(403) * density(25) 
  pd(26,25) = pd(26,25) - rrt(403) * density(26) 
  pd(26,26) = pd(26,26) - rrt(403) * density(25) 
  pd(37,25) = pd(37,25) + rrt(403) * density(26) 
  pd(37,26) = pd(37,26) + rrt(403) * density(25) 
  pd(02,26) = pd(02,26) + rrt(404) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(404) * density(26) * 2.0d0
  pd(25,26) = pd(25,26) + rrt(404) * density(34) 
  pd(25,34) = pd(25,34) + rrt(404) * density(26) 
  pd(26,26) = pd(26,26) - rrt(404) * density(34) 
  pd(26,34) = pd(26,34) - rrt(404) * density(26) 
  pd(34,26) = pd(34,26) - rrt(404) * density(34) 
  pd(34,34) = pd(34,34) - rrt(404) * density(26) 
  pd(02,26) = pd(02,26) + rrt(405) * density(35) * 2.0d0
  pd(02,35) = pd(02,35) + rrt(405) * density(26) * 2.0d0
  pd(15,26) = pd(15,26) + rrt(405) * density(35) 
  pd(15,35) = pd(15,35) + rrt(405) * density(26) 
  pd(26,26) = pd(26,26) - rrt(405) * density(35) 
  pd(26,35) = pd(26,35) - rrt(405) * density(26) 
  pd(35,26) = pd(35,26) - rrt(405) * density(35) 
  pd(35,35) = pd(35,35) - rrt(405) * density(26) 
  pd(09,26) = pd(09,26) + rrt(406) * density(40) 
  pd(09,40) = pd(09,40) + rrt(406) * density(26) 
  pd(15,26) = pd(15,26) + rrt(406) * density(40) 
  pd(15,40) = pd(15,40) + rrt(406) * density(26) 
  pd(25,26) = pd(25,26) + rrt(406) * density(40) 
  pd(25,40) = pd(25,40) + rrt(406) * density(26) 
  pd(26,26) = pd(26,26) - rrt(406) * density(40) 
  pd(26,40) = pd(26,40) - rrt(406) * density(26) 
  pd(40,26) = pd(40,26) - rrt(406) * density(40) 
  pd(40,40) = pd(40,40) - rrt(406) * density(26) 
  pd(09,25) = pd(09,25) + rrt(407) * density(28) 
  pd(09,28) = pd(09,28) + rrt(407) * density(25) 
  pd(25,25) = pd(25,25) - rrt(407) * density(28) 
  pd(25,28) = pd(25,28) - rrt(407) * density(25) 
  pd(28,25) = pd(28,25) - rrt(407) * density(28) 
  pd(28,28) = pd(28,28) - rrt(407) * density(25) 
  pd(37,25) = pd(37,25) + rrt(407) * density(28) 
  pd(37,28) = pd(37,28) + rrt(407) * density(25) 
  pd(09,29) = pd(09,29) + rrt(408) 
  pd(28,29) = pd(28,29) + rrt(408) 
  pd(29,29) = pd(29,29) - rrt(408) 
  pd(02,29) = pd(02,29) + rrt(409) * density(34) 
  pd(02,34) = pd(02,34) + rrt(409) * density(29) 
  pd(09,29) = pd(09,29) + rrt(409) * density(34) * 2.0d0
  pd(09,34) = pd(09,34) + rrt(409) * density(29) * 2.0d0
  pd(25,29) = pd(25,29) + rrt(409) * density(34) 
  pd(25,34) = pd(25,34) + rrt(409) * density(29) 
  pd(29,29) = pd(29,29) - rrt(409) * density(34) 
  pd(29,34) = pd(29,34) - rrt(409) * density(29) 
  pd(34,29) = pd(34,29) - rrt(409) * density(34) 
  pd(34,34) = pd(34,34) - rrt(409) * density(29) 
  pd(02,29) = pd(02,29) + rrt(410) * density(35) 
  pd(02,35) = pd(02,35) + rrt(410) * density(29) 
  pd(09,29) = pd(09,29) + rrt(410) * density(35) * 2.0d0
  pd(09,35) = pd(09,35) + rrt(410) * density(29) * 2.0d0
  pd(15,29) = pd(15,29) + rrt(410) * density(35) 
  pd(15,35) = pd(15,35) + rrt(410) * density(29) 
  pd(29,29) = pd(29,29) - rrt(410) * density(35) 
  pd(29,35) = pd(29,35) - rrt(410) * density(29) 
  pd(35,29) = pd(35,29) - rrt(410) * density(35) 
  pd(35,35) = pd(35,35) - rrt(410) * density(29) 
  pd(09,29) = pd(09,29) + rrt(411) * density(40) * 2.0d0
  pd(09,40) = pd(09,40) + rrt(411) * density(29) * 2.0d0
  pd(15,29) = pd(15,29) + rrt(411) * density(40) 
  pd(15,40) = pd(15,40) + rrt(411) * density(29) 
  pd(29,29) = pd(29,29) - rrt(411) * density(40) 
  pd(29,40) = pd(29,40) - rrt(411) * density(29) 
  pd(40,29) = pd(40,29) - rrt(411) * density(40) 
  pd(40,40) = pd(40,40) - rrt(411) * density(29) 
  pd(02,30) = pd(02,30) + rrt(412) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(412) * density(30) * 2.0d0
  pd(09,30) = pd(09,30) + rrt(412) * density(34) 
  pd(09,34) = pd(09,34) + rrt(412) * density(30) 
  pd(25,30) = pd(25,30) + rrt(412) * density(34) 
  pd(25,34) = pd(25,34) + rrt(412) * density(30) 
  pd(30,30) = pd(30,30) - rrt(412) * density(34) 
  pd(30,34) = pd(30,34) - rrt(412) * density(30) 
  pd(34,30) = pd(34,30) - rrt(412) * density(34) 
  pd(34,34) = pd(34,34) - rrt(412) * density(30) 
  pd(02,30) = pd(02,30) + rrt(413) * density(35) * 2.0d0
  pd(02,35) = pd(02,35) + rrt(413) * density(30) * 2.0d0
  pd(09,30) = pd(09,30) + rrt(413) * density(35) 
  pd(09,35) = pd(09,35) + rrt(413) * density(30) 
  pd(15,30) = pd(15,30) + rrt(413) * density(35) 
  pd(15,35) = pd(15,35) + rrt(413) * density(30) 
  pd(30,30) = pd(30,30) - rrt(413) * density(35) 
  pd(30,35) = pd(30,35) - rrt(413) * density(30) 
  pd(35,30) = pd(35,30) - rrt(413) * density(35) 
  pd(35,35) = pd(35,35) - rrt(413) * density(30) 
  pd(02,30) = pd(02,30) + rrt(414) * density(40) 
  pd(02,40) = pd(02,40) + rrt(414) * density(30) 
  pd(09,30) = pd(09,30) + rrt(414) * density(40) 
  pd(09,40) = pd(09,40) + rrt(414) * density(30) 
  pd(15,30) = pd(15,30) + rrt(414) * density(40) 
  pd(15,40) = pd(15,40) + rrt(414) * density(30) 
  pd(30,30) = pd(30,30) - rrt(414) * density(40) 
  pd(30,40) = pd(30,40) - rrt(414) * density(30) 
  pd(40,30) = pd(40,30) - rrt(414) * density(40) 
  pd(40,40) = pd(40,40) - rrt(414) * density(30) 
  pd(02,31) = pd(02,31) + rrt(415) 
  pd(26,31) = pd(26,31) + rrt(415) 
  pd(31,31) = pd(31,31) - rrt(415) 
  pd(02,31) = pd(02,31) + rrt(416) * density(34) * 3.0d0
  pd(02,34) = pd(02,34) + rrt(416) * density(31) * 3.0d0
  pd(25,31) = pd(25,31) + rrt(416) * density(34) 
  pd(25,34) = pd(25,34) + rrt(416) * density(31) 
  pd(31,31) = pd(31,31) - rrt(416) * density(34) 
  pd(31,34) = pd(31,34) - rrt(416) * density(31) 
  pd(34,31) = pd(34,31) - rrt(416) * density(34) 
  pd(34,34) = pd(34,34) - rrt(416) * density(31) 
  pd(02,31) = pd(02,31) + rrt(417) * density(35) * 3.0d0
  pd(02,35) = pd(02,35) + rrt(417) * density(31) * 3.0d0
  pd(15,31) = pd(15,31) + rrt(417) * density(35) 
  pd(15,35) = pd(15,35) + rrt(417) * density(31) 
  pd(31,31) = pd(31,31) - rrt(417) * density(35) 
  pd(31,35) = pd(31,35) - rrt(417) * density(31) 
  pd(35,31) = pd(35,31) - rrt(417) * density(35) 
  pd(35,35) = pd(35,35) - rrt(417) * density(31) 
  pd(02,31) = pd(02,31) + rrt(418) * density(40) * 2.0d0
  pd(02,40) = pd(02,40) + rrt(418) * density(31) * 2.0d0
  pd(15,31) = pd(15,31) + rrt(418) * density(40) 
  pd(15,40) = pd(15,40) + rrt(418) * density(31) 
  pd(31,31) = pd(31,31) - rrt(418) * density(40) 
  pd(31,40) = pd(31,40) - rrt(418) * density(31) 
  pd(40,31) = pd(40,31) - rrt(418) * density(40) 
  pd(40,40) = pd(40,40) - rrt(418) * density(31) 
  pd(02,34) = pd(02,34) + rrt(419) * density(36) 
  pd(02,36) = pd(02,36) + rrt(419) * density(34) 
  pd(15,34) = pd(15,34) + rrt(419) * density(36) 
  pd(15,36) = pd(15,36) + rrt(419) * density(34) 
  pd(25,34) = pd(25,34) + rrt(419) * density(36) 
  pd(25,36) = pd(25,36) + rrt(419) * density(34) 
  pd(34,34) = pd(34,34) - rrt(419) * density(36) 
  pd(34,36) = pd(34,36) - rrt(419) * density(34) 
  pd(36,34) = pd(36,34) - rrt(419) * density(36) 
  pd(36,36) = pd(36,36) - rrt(419) * density(34) 
  pd(02,25) = pd(02,25) + rrt(420) * density(34) 
  pd(02,34) = pd(02,34) + rrt(420) * density(25) 
  pd(25,25) = pd(25,25) - rrt(420) * density(34) 
  pd(25,34) = pd(25,34) - rrt(420) * density(25) 
  pd(34,25) = pd(34,25) - rrt(420) * density(34) 
  pd(34,34) = pd(34,34) - rrt(420) * density(25) 
  pd(40,25) = pd(40,25) + rrt(420) * density(34) 
  pd(40,34) = pd(40,34) + rrt(420) * density(25) 
  pd(02,35) = pd(02,35) + rrt(421) * density(36) 
  pd(02,36) = pd(02,36) + rrt(421) * density(35) 
  pd(15,35) = pd(15,35) + rrt(421) * density(36) * 2.0d0
  pd(15,36) = pd(15,36) + rrt(421) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(421) * density(36) 
  pd(35,36) = pd(35,36) - rrt(421) * density(35) 
  pd(36,35) = pd(36,35) - rrt(421) * density(36) 
  pd(36,36) = pd(36,36) - rrt(421) * density(35) 
  pd(15,25) = pd(15,25) + rrt(422) * density(35) 
  pd(15,35) = pd(15,35) + rrt(422) * density(25) 
  pd(25,25) = pd(25,25) - rrt(422) * density(35) 
  pd(25,35) = pd(25,35) - rrt(422) * density(25) 
  pd(34,25) = pd(34,25) + rrt(422) * density(35) 
  pd(34,35) = pd(34,35) + rrt(422) * density(25) 
  pd(35,25) = pd(35,25) - rrt(422) * density(35) 
  pd(35,35) = pd(35,35) - rrt(422) * density(25) 
  pd(02,25) = pd(02,25) + rrt(423) * density(35) 
  pd(02,35) = pd(02,35) + rrt(423) * density(25) 
  pd(15,25) = pd(15,25) + rrt(423) * density(35) 
  pd(15,35) = pd(15,35) + rrt(423) * density(25) 
  pd(25,25) = pd(25,25) - rrt(423) * density(35) 
  pd(25,35) = pd(25,35) - rrt(423) * density(25) 
  pd(35,25) = pd(35,25) - rrt(423) * density(35) 
  pd(35,35) = pd(35,35) - rrt(423) * density(25) 
  pd(39,25) = pd(39,25) + rrt(423) * density(35) 
  pd(39,35) = pd(39,35) + rrt(423) * density(25) 
  pd(02,25) = pd(02,25) + rrt(424) * density(35) 
  pd(02,35) = pd(02,35) + rrt(424) * density(25) 
  pd(25,25) = pd(25,25) - rrt(424) * density(35) 
  pd(25,35) = pd(25,35) - rrt(424) * density(25) 
  pd(35,25) = pd(35,25) - rrt(424) * density(35) 
  pd(35,35) = pd(35,35) - rrt(424) * density(25) 
  pd(41,25) = pd(41,25) + rrt(424) * density(35) 
  pd(41,35) = pd(41,35) + rrt(424) * density(25) 
  pd(02,24) = pd(02,24) + rrt(425) * density(35) 
  pd(02,35) = pd(02,35) + rrt(425) * density(24) 
  pd(15,24) = pd(15,24) + rrt(425) * density(35) 
  pd(15,35) = pd(15,35) + rrt(425) * density(24) 
  pd(24,24) = pd(24,24) - rrt(425) * density(35) 
  pd(24,35) = pd(24,35) - rrt(425) * density(24) 
  pd(35,24) = pd(35,24) - rrt(425) * density(35) 
  pd(35,35) = pd(35,35) - rrt(425) * density(24) 
  pd(41,24) = pd(41,24) + rrt(425) * density(35) 
  pd(41,35) = pd(41,35) + rrt(425) * density(24) 
  pd(24,24) = pd(24,24) - rrt(426) * density(39) 
  pd(24,39) = pd(24,39) - rrt(426) * density(24) 
  pd(25,24) = pd(25,24) + rrt(426) * density(39) 
  pd(25,39) = pd(25,39) + rrt(426) * density(24) 
  pd(39,24) = pd(39,24) - rrt(426) * density(39) 
  pd(39,39) = pd(39,39) - rrt(426) * density(24) 
  pd(41,24) = pd(41,24) + rrt(426) * density(39) 
  pd(41,39) = pd(41,39) + rrt(426) * density(24) 
  pd(01,24) = pd(01,24) + rrt(427) * density(39) 
  pd(01,39) = pd(01,39) + rrt(427) * density(24) 
  pd(15,24) = pd(15,24) + rrt(427) * density(39) * 2.0d0
  pd(15,39) = pd(15,39) + rrt(427) * density(24) * 2.0d0
  pd(24,24) = pd(24,24) - rrt(427) * density(39) 
  pd(24,39) = pd(24,39) - rrt(427) * density(24) 
  pd(39,24) = pd(39,24) - rrt(427) * density(39) 
  pd(39,39) = pd(39,39) - rrt(427) * density(24) 
  pd(15,24) = pd(15,24) + rrt(428) * density(40) 
  pd(15,40) = pd(15,40) + rrt(428) * density(24) 
  pd(24,24) = pd(24,24) - rrt(428) * density(40) 
  pd(24,40) = pd(24,40) - rrt(428) * density(24) 
  pd(40,24) = pd(40,24) - rrt(428) * density(40) 
  pd(40,40) = pd(40,40) - rrt(428) * density(24) 
  pd(41,24) = pd(41,24) + rrt(428) * density(40) 
  pd(41,40) = pd(41,40) + rrt(428) * density(24) 
  pd(01,24) = pd(01,24) + rrt(429) * density(41) 
  pd(01,41) = pd(01,41) + rrt(429) * density(24) 
  pd(15,24) = pd(15,24) + rrt(429) * density(41) * 3.0d0
  pd(15,41) = pd(15,41) + rrt(429) * density(24) * 3.0d0
  pd(24,24) = pd(24,24) - rrt(429) * density(41) 
  pd(24,41) = pd(24,41) - rrt(429) * density(24) 
  pd(41,24) = pd(41,24) - rrt(429) * density(41) 
  pd(41,41) = pd(41,41) - rrt(429) * density(24) 
  pd(15,24) = pd(15,24) + rrt(430) * density(37) 
  pd(15,37) = pd(15,37) + rrt(430) * density(24) 
  pd(24,24) = pd(24,24) - rrt(430) * density(37) 
  pd(24,37) = pd(24,37) - rrt(430) * density(24) 
  pd(36,24) = pd(36,24) + rrt(430) * density(37) 
  pd(36,37) = pd(36,37) + rrt(430) * density(24) 
  pd(37,24) = pd(37,24) - rrt(430) * density(37) 
  pd(37,37) = pd(37,37) - rrt(430) * density(24) 
  pd(25,25) = pd(25,25) - rrt(431) * density(37) 
  pd(25,37) = pd(25,37) - rrt(431) * density(25) 
  pd(36,25) = pd(36,25) + rrt(431) * density(37) 
  pd(36,37) = pd(36,37) + rrt(431) * density(25) 
  pd(37,25) = pd(37,25) - rrt(431) * density(37) 
  pd(37,37) = pd(37,37) - rrt(431) * density(25) 
  pd(01,25) = pd(01,25) + rrt(432) * density(39) 
  pd(01,39) = pd(01,39) + rrt(432) * density(25) 
  pd(15,25) = pd(15,25) + rrt(432) * density(39) 
  pd(15,39) = pd(15,39) + rrt(432) * density(25) 
  pd(25,25) = pd(25,25) - rrt(432) * density(39) 
  pd(25,39) = pd(25,39) - rrt(432) * density(25) 
  pd(39,25) = pd(39,25) - rrt(432) * density(39) 
  pd(39,39) = pd(39,39) - rrt(432) * density(25) 
  pd(15,25) = pd(15,25) + rrt(433) * density(40) 
  pd(15,40) = pd(15,40) + rrt(433) * density(25) 
  pd(25,25) = pd(25,25) - rrt(433) * density(40) 
  pd(25,40) = pd(25,40) - rrt(433) * density(25) 
  pd(39,25) = pd(39,25) + rrt(433) * density(40) 
  pd(39,40) = pd(39,40) + rrt(433) * density(25) 
  pd(40,25) = pd(40,25) - rrt(433) * density(40) 
  pd(40,40) = pd(40,40) - rrt(433) * density(25) 
  pd(01,25) = pd(01,25) + rrt(434) * density(40) 
  pd(01,40) = pd(01,40) + rrt(434) * density(25) 
  pd(24,25) = pd(24,25) + rrt(434) * density(40) 
  pd(24,40) = pd(24,40) + rrt(434) * density(25) 
  pd(25,25) = pd(25,25) - rrt(434) * density(40) 
  pd(25,40) = pd(25,40) - rrt(434) * density(25) 
  pd(40,25) = pd(40,25) - rrt(434) * density(40) 
  pd(40,40) = pd(40,40) - rrt(434) * density(25) 
  pd(24,25) = pd(24,25) + rrt(435) * density(41) 
  pd(24,41) = pd(24,41) + rrt(435) * density(25) 
  pd(25,25) = pd(25,25) - rrt(435) * density(41) 
  pd(25,41) = pd(25,41) - rrt(435) * density(25) 
  pd(39,25) = pd(39,25) + rrt(435) * density(41) 
  pd(39,41) = pd(39,41) + rrt(435) * density(25) 
  pd(41,25) = pd(41,25) - rrt(435) * density(41) 
  pd(41,41) = pd(41,41) - rrt(435) * density(25) 
  pd(01,25) = pd(01,25) + rrt(436) * density(41) 
  pd(01,41) = pd(01,41) + rrt(436) * density(25) 
  pd(15,25) = pd(15,25) + rrt(436) * density(41) * 2.0d0
  pd(15,41) = pd(15,41) + rrt(436) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(436) * density(41) 
  pd(25,41) = pd(25,41) - rrt(436) * density(25) 
  pd(41,25) = pd(41,25) - rrt(436) * density(41) 
  pd(41,41) = pd(41,41) - rrt(436) * density(25) 
  pd(15,25) = pd(15,25) + rrt(437) * density(41) 
  pd(15,41) = pd(15,41) + rrt(437) * density(25) 
  pd(25,25) = pd(25,25) - rrt(437) * density(41) 
  pd(25,41) = pd(25,41) - rrt(437) * density(25) 
  pd(40,25) = pd(40,25) + rrt(437) * density(41) 
  pd(40,41) = pd(40,41) + rrt(437) * density(25) 
  pd(41,25) = pd(41,25) - rrt(437) * density(41) 
  pd(41,41) = pd(41,41) - rrt(437) * density(25) 
  pd(15,25) = pd(15,25) + rrt(438) * density(42) 
  pd(15,42) = pd(15,42) + rrt(438) * density(25) 
  pd(25,25) = pd(25,25) - rrt(438) * density(42) 
  pd(25,42) = pd(25,42) - rrt(438) * density(25) 
  pd(41,25) = pd(41,25) + rrt(438) * density(42) 
  pd(41,42) = pd(41,42) + rrt(438) * density(25) 
  pd(42,25) = pd(42,25) - rrt(438) * density(42) 
  pd(42,42) = pd(42,42) - rrt(438) * density(25) 
  pd(15,25) = pd(15,25) + rrt(439) * density(42) * 2.0d0
  pd(15,42) = pd(15,42) + rrt(439) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(439) * density(42) 
  pd(25,42) = pd(25,42) - rrt(439) * density(25) 
  pd(39,25) = pd(39,25) + rrt(439) * density(42) 
  pd(39,42) = pd(39,42) + rrt(439) * density(25) 
  pd(42,25) = pd(42,25) - rrt(439) * density(42) 
  pd(42,42) = pd(42,42) - rrt(439) * density(25) 
  pd(24,25) = pd(24,25) + rrt(440) * density(42) 
  pd(24,42) = pd(24,42) + rrt(440) * density(25) 
  pd(25,25) = pd(25,25) - rrt(440) * density(42) 
  pd(25,42) = pd(25,42) - rrt(440) * density(25) 
  pd(40,25) = pd(40,25) + rrt(440) * density(42) 
  pd(40,42) = pd(40,42) + rrt(440) * density(25) 
  pd(42,25) = pd(42,25) - rrt(440) * density(42) 
  pd(42,42) = pd(42,42) - rrt(440) * density(25) 
  pd(24,37) = pd(24,37) + rrt(441) * density(40) 
  pd(24,40) = pd(24,40) + rrt(441) * density(37) 
  pd(37,37) = pd(37,37) - rrt(441) * density(40) 
  pd(37,40) = pd(37,40) - rrt(441) * density(37) 
  pd(40,37) = pd(40,37) - rrt(441) * density(40) 
  pd(40,40) = pd(40,40) - rrt(441) * density(37) 
  pd(15,37) = pd(15,37) + rrt(442) * density(40) 
  pd(15,40) = pd(15,40) + rrt(442) * density(37) 
  pd(25,37) = pd(25,37) + rrt(442) * density(40) 
  pd(25,40) = pd(25,40) + rrt(442) * density(37) 
  pd(37,37) = pd(37,37) - rrt(442) * density(40) 
  pd(37,40) = pd(37,40) - rrt(442) * density(37) 
  pd(40,37) = pd(40,37) - rrt(442) * density(40) 
  pd(40,40) = pd(40,40) - rrt(442) * density(37) 
  pd(15,36) = pd(15,36) + rrt(443) * density(40) * 2.0d0
  pd(15,40) = pd(15,40) + rrt(443) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(443) * density(40) 
  pd(36,40) = pd(36,40) - rrt(443) * density(36) 
  pd(40,36) = pd(40,36) - rrt(443) * density(40) 
  pd(40,40) = pd(40,40) - rrt(443) * density(36) 
  pd(15,36) = pd(15,36) + rrt(444) * density(40) 
  pd(15,40) = pd(15,40) + rrt(444) * density(36) 
  pd(25,36) = pd(25,36) + rrt(444) * density(40) * 2.0d0
  pd(25,40) = pd(25,40) + rrt(444) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(444) * density(40) 
  pd(36,40) = pd(36,40) - rrt(444) * density(36) 
  pd(40,36) = pd(40,36) - rrt(444) * density(40) 
  pd(40,40) = pd(40,40) - rrt(444) * density(36) 
  pd(15,36) = pd(15,36) + rrt(445) * density(40) * 2.0d0
  pd(15,40) = pd(15,40) + rrt(445) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(445) * density(40) 
  pd(36,40) = pd(36,40) - rrt(445) * density(36) 
  pd(40,36) = pd(40,36) - rrt(445) * density(40) 
  pd(40,40) = pd(40,40) - rrt(445) * density(36) 
  pd(01,40) = pd(01,40) + rrt(446) 
  pd(15,40) = pd(15,40) + rrt(446) 
  pd(40,40) = pd(40,40) - rrt(446) 
  pd(15,36) = pd(15,36) + rrt(447) * density(41) 
  pd(15,41) = pd(15,41) + rrt(447) * density(36) 
  pd(24,36) = pd(24,36) + rrt(447) * density(41) 
  pd(24,41) = pd(24,41) + rrt(447) * density(36) 
  pd(36,36) = pd(36,36) - rrt(447) * density(41) 
  pd(36,41) = pd(36,41) - rrt(447) * density(36) 
  pd(41,36) = pd(41,36) - rrt(447) * density(41) 
  pd(41,41) = pd(41,41) - rrt(447) * density(36) 
  pd(24,36) = pd(24,36) + rrt(448) * density(41) 
  pd(24,41) = pd(24,41) + rrt(448) * density(36) 
  pd(25,36) = pd(25,36) + rrt(448) * density(41) * 2.0d0
  pd(25,41) = pd(25,41) + rrt(448) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(448) * density(41) 
  pd(36,41) = pd(36,41) - rrt(448) * density(36) 
  pd(41,36) = pd(41,36) - rrt(448) * density(41) 
  pd(41,41) = pd(41,41) - rrt(448) * density(36) 
  pd(24,37) = pd(24,37) + rrt(449) * density(41) 
  pd(24,41) = pd(24,41) + rrt(449) * density(37) 
  pd(25,37) = pd(25,37) + rrt(449) * density(41) 
  pd(25,41) = pd(25,41) + rrt(449) * density(37) 
  pd(37,37) = pd(37,37) - rrt(449) * density(41) 
  pd(37,41) = pd(37,41) - rrt(449) * density(37) 
  pd(41,37) = pd(41,37) - rrt(449) * density(41) 
  pd(41,41) = pd(41,41) - rrt(449) * density(37) 
  pd(01,41) = pd(01,41) + rrt(450) 
  pd(24,41) = pd(24,41) + rrt(450) 
  pd(41,41) = pd(41,41) - rrt(450) 
  pd(25,37) = pd(25,37) + rrt(451) * density(39) * 2.0d0
  pd(25,39) = pd(25,39) + rrt(451) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(451) * density(39) 
  pd(37,39) = pd(37,39) - rrt(451) * density(37) 
  pd(39,37) = pd(39,37) - rrt(451) * density(39) 
  pd(39,39) = pd(39,39) - rrt(451) * density(37) 
  pd(15,37) = pd(15,37) + rrt(452) * density(39) 
  pd(15,39) = pd(15,39) + rrt(452) * density(37) 
  pd(37,37) = pd(37,37) - rrt(452) * density(39) 
  pd(37,39) = pd(37,39) - rrt(452) * density(37) 
  pd(39,37) = pd(39,37) - rrt(452) * density(39) 
  pd(39,39) = pd(39,39) - rrt(452) * density(37) 
  pd(15,36) = pd(15,36) + rrt(453) * density(39) 
  pd(15,39) = pd(15,39) + rrt(453) * density(36) 
  pd(25,36) = pd(25,36) + rrt(453) * density(39) 
  pd(25,39) = pd(25,39) + rrt(453) * density(36) 
  pd(36,36) = pd(36,36) - rrt(453) * density(39) 
  pd(36,39) = pd(36,39) - rrt(453) * density(36) 
  pd(39,36) = pd(39,36) - rrt(453) * density(39) 
  pd(39,39) = pd(39,39) - rrt(453) * density(36) 
  pd(25,36) = pd(25,36) + rrt(454) * density(39) * 3.0d0
  pd(25,39) = pd(25,39) + rrt(454) * density(36) * 3.0d0
  pd(36,36) = pd(36,36) - rrt(454) * density(39) 
  pd(36,39) = pd(36,39) - rrt(454) * density(36) 
  pd(39,36) = pd(39,36) - rrt(454) * density(39) 
  pd(39,39) = pd(39,39) - rrt(454) * density(36) 
  pd(24,36) = pd(24,36) + rrt(455) * density(39) 
  pd(24,39) = pd(24,39) + rrt(455) * density(36) 
  pd(36,36) = pd(36,36) - rrt(455) * density(39) 
  pd(36,39) = pd(36,39) - rrt(455) * density(36) 
  pd(39,36) = pd(39,36) - rrt(455) * density(39) 
  pd(39,39) = pd(39,39) - rrt(455) * density(36) 
  pd(01,39) = pd(01,39) + rrt(456) 
  pd(25,39) = pd(25,39) + rrt(456) 
  pd(39,39) = pd(39,39) - rrt(456) 
  pd(15,42) = pd(15,42) + rrt(457) 
  pd(40,42) = pd(40,42) + rrt(457) 
  pd(42,42) = pd(42,42) - rrt(457) 
  pd(15,38) = pd(15,38) + rrt(458) 
  pd(36,38) = pd(36,38) + rrt(458) 
  pd(38,38) = pd(38,38) - rrt(458) 
  pd(02,02) = pd(02,02) + rrt(459) * density(03) 
  pd(02,03) = pd(02,03) + rrt(459) * density(02) 
  pd(03,02) = pd(03,02) - rrt(459) * density(03) 
  pd(03,03) = pd(03,03) - rrt(459) * density(02) 
  pd(02,03) = pd(02,03) + rrt(460) * density(09) 
  pd(02,09) = pd(02,09) + rrt(460) * density(03) 
  pd(03,03) = pd(03,03) - rrt(460) * density(09) 
  pd(03,09) = pd(03,09) - rrt(460) * density(03) 
  pd(02,03) = pd(02,03) + rrt(461) * density(15) 
  pd(02,15) = pd(02,15) + rrt(461) * density(03) 
  pd(03,03) = pd(03,03) - rrt(461) * density(15) 
  pd(03,15) = pd(03,15) - rrt(461) * density(03) 
  pd(02,02) = pd(02,02) + rrt(462) * density(04) 
  pd(02,04) = pd(02,04) + rrt(462) * density(02) 
  pd(04,02) = pd(04,02) - rrt(462) * density(04) 
  pd(04,04) = pd(04,04) - rrt(462) * density(02) 
  pd(02,04) = pd(02,04) + rrt(463) * density(09) 
  pd(02,09) = pd(02,09) + rrt(463) * density(04) 
  pd(04,04) = pd(04,04) - rrt(463) * density(09) 
  pd(04,09) = pd(04,09) - rrt(463) * density(04) 
  pd(02,04) = pd(02,04) + rrt(464) * density(15) 
  pd(02,15) = pd(02,15) + rrt(464) * density(04) 
  pd(04,04) = pd(04,04) - rrt(464) * density(15) 
  pd(04,15) = pd(04,15) - rrt(464) * density(04) 
  pd(03,02) = pd(03,02) + rrt(465) * density(04) 
  pd(03,04) = pd(03,04) + rrt(465) * density(02) 
  pd(04,02) = pd(04,02) - rrt(465) * density(04) 
  pd(04,04) = pd(04,04) - rrt(465) * density(02) 
  pd(03,04) = pd(03,04) + rrt(466) * density(09) 
  pd(03,09) = pd(03,09) + rrt(466) * density(04) 
  pd(04,04) = pd(04,04) - rrt(466) * density(09) 
  pd(04,09) = pd(04,09) - rrt(466) * density(04) 
  pd(03,04) = pd(03,04) + rrt(467) * density(15) 
  pd(03,15) = pd(03,15) + rrt(467) * density(04) 
  pd(04,04) = pd(04,04) - rrt(467) * density(15) 
  pd(04,15) = pd(04,15) - rrt(467) * density(04) 
  pd(04,02) = pd(04,02) + rrt(468) * density(05) 
  pd(04,05) = pd(04,05) + rrt(468) * density(02) 
  pd(05,02) = pd(05,02) - rrt(468) * density(05) 
  pd(05,05) = pd(05,05) - rrt(468) * density(02) 
  pd(04,05) = pd(04,05) + rrt(469) * density(09) 
  pd(04,09) = pd(04,09) + rrt(469) * density(05) 
  pd(05,05) = pd(05,05) - rrt(469) * density(09) 
  pd(05,09) = pd(05,09) - rrt(469) * density(05) 
  pd(04,05) = pd(04,05) + rrt(470) * density(15) 
  pd(04,15) = pd(04,15) + rrt(470) * density(05) 
  pd(05,05) = pd(05,05) - rrt(470) * density(15) 
  pd(05,15) = pd(05,15) - rrt(470) * density(05) 
  pd(05,02) = pd(05,02) - rrt(471) * density(05) 
  pd(05,05) = pd(05,05) - rrt(471) * density(02) 
  pd(06,02) = pd(06,02) + rrt(471) * density(05) 
  pd(06,05) = pd(06,05) + rrt(471) * density(02) 
  pd(05,05) = pd(05,05) - rrt(472) * density(09) 
  pd(05,09) = pd(05,09) - rrt(472) * density(05) 
  pd(06,05) = pd(06,05) + rrt(472) * density(09) 
  pd(06,09) = pd(06,09) + rrt(472) * density(05) 
  pd(05,05) = pd(05,05) - rrt(473) * density(15) 
  pd(05,15) = pd(05,15) - rrt(473) * density(05) 
  pd(06,05) = pd(06,05) + rrt(473) * density(15) 
  pd(06,15) = pd(06,15) + rrt(473) * density(05) 
  pd(02,02) = pd(02,02) - rrt(474) * density(05) 
  pd(02,05) = pd(02,05) - rrt(474) * density(02) 
  pd(03,02) = pd(03,02) + rrt(474) * density(05) 
  pd(03,05) = pd(03,05) + rrt(474) * density(02) 
  pd(04,02) = pd(04,02) + rrt(474) * density(05) 
  pd(04,05) = pd(04,05) + rrt(474) * density(02) 
  pd(05,02) = pd(05,02) - rrt(474) * density(05) 
  pd(05,05) = pd(05,05) - rrt(474) * density(02) 
  pd(03,02) = pd(03,02) + rrt(475) * density(05) 
  pd(03,05) = pd(03,05) + rrt(475) * density(02) 
  pd(05,02) = pd(05,02) - rrt(475) * density(05) 
  pd(05,05) = pd(05,05) - rrt(475) * density(02) 
  pd(03,05) = pd(03,05) + rrt(476) * density(09) 
  pd(03,09) = pd(03,09) + rrt(476) * density(05) 
  pd(05,05) = pd(05,05) - rrt(476) * density(09) 
  pd(05,09) = pd(05,09) - rrt(476) * density(05) 
  pd(03,05) = pd(03,05) + rrt(477) * density(15) 
  pd(03,15) = pd(03,15) + rrt(477) * density(05) 
  pd(05,05) = pd(05,05) - rrt(477) * density(15) 
  pd(05,15) = pd(05,15) - rrt(477) * density(05) 
  pd(04,02) = pd(04,02) + rrt(478) * density(06) 
  pd(04,06) = pd(04,06) + rrt(478) * density(02) 
  pd(06,02) = pd(06,02) - rrt(478) * density(06) 
  pd(06,06) = pd(06,06) - rrt(478) * density(02) 
  pd(04,06) = pd(04,06) + rrt(479) * density(09) 
  pd(04,09) = pd(04,09) + rrt(479) * density(06) 
  pd(06,06) = pd(06,06) - rrt(479) * density(09) 
  pd(06,09) = pd(06,09) - rrt(479) * density(06) 
  pd(04,06) = pd(04,06) + rrt(480) * density(15) 
  pd(04,15) = pd(04,15) + rrt(480) * density(06) 
  pd(06,06) = pd(06,06) - rrt(480) * density(15) 
  pd(06,15) = pd(06,15) - rrt(480) * density(06) 
  pd(03,02) = pd(03,02) + rrt(481) * density(06) 
  pd(03,06) = pd(03,06) + rrt(481) * density(02) 
  pd(06,02) = pd(06,02) - rrt(481) * density(06) 
  pd(06,06) = pd(06,06) - rrt(481) * density(02) 
  pd(03,06) = pd(03,06) + rrt(482) * density(09) 
  pd(03,09) = pd(03,09) + rrt(482) * density(06) 
  pd(06,06) = pd(06,06) - rrt(482) * density(09) 
  pd(06,09) = pd(06,09) - rrt(482) * density(06) 
  pd(03,06) = pd(03,06) + rrt(483) * density(15) 
  pd(03,15) = pd(03,15) + rrt(483) * density(06) 
  pd(06,06) = pd(06,06) - rrt(483) * density(15) 
  pd(06,15) = pd(06,15) - rrt(483) * density(06) 
  pd(09,02) = pd(09,02) + rrt(484) * density(10) 
  pd(09,10) = pd(09,10) + rrt(484) * density(02) 
  pd(10,02) = pd(10,02) - rrt(484) * density(10) 
  pd(10,10) = pd(10,10) - rrt(484) * density(02) 
  pd(09,09) = pd(09,09) + rrt(485) * density(10) 
  pd(09,10) = pd(09,10) + rrt(485) * density(09) 
  pd(10,09) = pd(10,09) - rrt(485) * density(10) 
  pd(10,10) = pd(10,10) - rrt(485) * density(09) 
  pd(09,10) = pd(09,10) + rrt(486) * density(15) 
  pd(09,15) = pd(09,15) + rrt(486) * density(10) 
  pd(10,10) = pd(10,10) - rrt(486) * density(15) 
  pd(10,15) = pd(10,15) - rrt(486) * density(10) 
  pd(15,02) = pd(15,02) + rrt(487) * density(16) 
  pd(15,16) = pd(15,16) + rrt(487) * density(02) 
  pd(16,02) = pd(16,02) - rrt(487) * density(16) 
  pd(16,16) = pd(16,16) - rrt(487) * density(02) 
  pd(15,09) = pd(15,09) + rrt(488) * density(16) 
  pd(15,16) = pd(15,16) + rrt(488) * density(09) 
  pd(16,09) = pd(16,09) - rrt(488) * density(16) 
  pd(16,16) = pd(16,16) - rrt(488) * density(09) 
  pd(15,15) = pd(15,15) + rrt(489) * density(16) 
  pd(15,16) = pd(15,16) + rrt(489) * density(15) 
  pd(16,15) = pd(16,15) - rrt(489) * density(16) 
  pd(16,16) = pd(16,16) - rrt(489) * density(15) 
  pd(15,02) = pd(15,02) + rrt(490) * density(17) 
  pd(15,17) = pd(15,17) + rrt(490) * density(02) 
  pd(17,02) = pd(17,02) - rrt(490) * density(17) 
  pd(17,17) = pd(17,17) - rrt(490) * density(02) 
  pd(15,09) = pd(15,09) + rrt(491) * density(17) 
  pd(15,17) = pd(15,17) + rrt(491) * density(09) 
  pd(17,09) = pd(17,09) - rrt(491) * density(17) 
  pd(17,17) = pd(17,17) - rrt(491) * density(09) 
  pd(15,15) = pd(15,15) + rrt(492) * density(17) 
  pd(15,17) = pd(15,17) + rrt(492) * density(15) 
  pd(17,15) = pd(17,15) - rrt(492) * density(17) 
  pd(17,17) = pd(17,17) - rrt(492) * density(15) 
  pd(15,02) = pd(15,02) + rrt(493) * density(18) 
  pd(15,18) = pd(15,18) + rrt(493) * density(02) 
  pd(18,02) = pd(18,02) - rrt(493) * density(18) 
  pd(18,18) = pd(18,18) - rrt(493) * density(02) 
  pd(15,09) = pd(15,09) + rrt(494) * density(18) 
  pd(15,18) = pd(15,18) + rrt(494) * density(09) 
  pd(18,09) = pd(18,09) - rrt(494) * density(18) 
  pd(18,18) = pd(18,18) - rrt(494) * density(09) 
  pd(15,15) = pd(15,15) + rrt(495) * density(18) 
  pd(15,18) = pd(15,18) + rrt(495) * density(15) 
  pd(18,15) = pd(18,15) - rrt(495) * density(18) 
  pd(18,18) = pd(18,18) - rrt(495) * density(15) 
  pd(02,07) = pd(02,07) + rrt(496) 
  pd(07,07) = pd(07,07) - rrt(496) 
  pd(02,08) = pd(02,08) + rrt(497) 
  pd(08,08) = pd(08,08) - rrt(497) 
  pd(09,11) = pd(09,11) + rrt(498) 
  pd(11,11) = pd(11,11) - rrt(498) 
  pd(09,12) = pd(09,12) + rrt(499) 
  pd(12,12) = pd(12,12) - rrt(499) 
  pd(09,13) = pd(09,13) + rrt(500) 
  pd(13,13) = pd(13,13) - rrt(500) 
  pd(09,14) = pd(09,14) + rrt(501) 
  pd(14,14) = pd(14,14) - rrt(501) 
  pd(15,19) = pd(15,19) + rrt(502) 
  pd(19,19) = pd(19,19) - rrt(502) 
  pd(15,20) = pd(15,20) + rrt(503) 
  pd(20,20) = pd(20,20) - rrt(503) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(43,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(43,:) = pd(43,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: Te
  double precision :: ANY_NEUTRAL
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  Te  = ZDPlasKin_cfg(4)
  call ZDPlasKin_get_density_total(ALL_NEUTRAL=ANY_NEUTRAL)
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = 0.05*0.65*1.4*1.0E-11*TE**(0.5)*(1+1.3*1.0E-5*TE)*EXP(-1.5*1.0E5/TE)
  rrt(006) = bolsig_rates(bolsig_pointer(5))
  rrt(007) = bolsig_rates(bolsig_pointer(6))
  rrt(008) = bolsig_rates(bolsig_pointer(7))
  rrt(009) = bolsig_rates(bolsig_pointer(8))
  rrt(010) = bolsig_rates(bolsig_pointer(9))
  rrt(011) = bolsig_rates(bolsig_pointer(10))
  rrt(012) = bolsig_rates(bolsig_pointer(11))
  rrt(013) = bolsig_rates(bolsig_pointer(12))
  rrt(014) = bolsig_rates(bolsig_pointer(13))
  rrt(015) = bolsig_rates(bolsig_pointer(14))
  rrt(016) = bolsig_rates(bolsig_pointer(15))
  rrt(017) = bolsig_rates(bolsig_pointer(16))
  rrt(018) = bolsig_rates(bolsig_pointer(17))
  rrt(019) = bolsig_rates(bolsig_pointer(18))
  rrt(020) = bolsig_rates(bolsig_pointer(19))
  rrt(021) = bolsig_rates(bolsig_pointer(20))
  rrt(022) = bolsig_rates(bolsig_pointer(21))
  rrt(023) = bolsig_rates(bolsig_pointer(22))
  rrt(024) = bolsig_rates(bolsig_pointer(23))
  rrt(025) = bolsig_rates(bolsig_pointer(24))
  rrt(026) = bolsig_rates(bolsig_pointer(25))
  rrt(027) = bolsig_rates(bolsig_pointer(26))
  rrt(028) = bolsig_rates(bolsig_pointer(27))
  rrt(029) = bolsig_rates(bolsig_pointer(28))
  rrt(030) = bolsig_rates(bolsig_pointer(29))
  rrt(031) = bolsig_rates(bolsig_pointer(30))
  rrt(032) = bolsig_rates(bolsig_pointer(31))
  rrt(033) = bolsig_rates(bolsig_pointer(32))
  rrt(034) = bolsig_rates(bolsig_pointer(33))
  rrt(035) = bolsig_rates(bolsig_pointer(34))
  rrt(036) = bolsig_rates(bolsig_pointer(35))
  rrt(037) = bolsig_rates(bolsig_pointer(36))
  rrt(038) = bolsig_rates(bolsig_pointer(37))
  rrt(039) = bolsig_rates(bolsig_pointer(38))
  rrt(040) = bolsig_rates(bolsig_pointer(39))
  rrt(041) = bolsig_rates(bolsig_pointer(40))
  rrt(042) = bolsig_rates(bolsig_pointer(41))
  rrt(043) = bolsig_rates(bolsig_pointer(42))
  rrt(044) = bolsig_rates(bolsig_pointer(43))
  rrt(045) = bolsig_rates(bolsig_pointer(44))
  rrt(046) = bolsig_rates(bolsig_pointer(45))
  rrt(047) = bolsig_rates(bolsig_pointer(46))
  rrt(048) = rrt(5)
  rrt(049) = rrt(5)
  rrt(050) = rrt(5)
  rrt(051) = rrt(5)
  rrt(052) = bolsig_rates(bolsig_pointer(47))
  rrt(053) = bolsig_rates(bolsig_pointer(48))
  rrt(054) = bolsig_rates(bolsig_pointer(49))
  rrt(055) = bolsig_rates(bolsig_pointer(50))
  rrt(056) = bolsig_rates(bolsig_pointer(51))
  rrt(057) = bolsig_rates(bolsig_pointer(52))
  rrt(058) = bolsig_rates(bolsig_pointer(53))
  rrt(059) = bolsig_rates(bolsig_pointer(54))
  rrt(060) = bolsig_rates(bolsig_pointer(55))
  rrt(061) = bolsig_rates(bolsig_pointer(56))
  rrt(062) = bolsig_rates(bolsig_pointer(57))
  rrt(063) = bolsig_rates(bolsig_pointer(58))
  rrt(064) = bolsig_rates(bolsig_pointer(59))
  rrt(065) = bolsig_rates(bolsig_pointer(60))
  rrt(066) = rrt(5)
  rrt(067) = rrt(5)
  rrt(068) = bolsig_rates(bolsig_pointer(61))
  rrt(069) = bolsig_rates(bolsig_pointer(62))
  rrt(070) = bolsig_rates(bolsig_pointer(63))
  rrt(071) = bolsig_rates(bolsig_pointer(64))
  rrt(072) = bolsig_rates(bolsig_pointer(65))
  rrt(073) = bolsig_rates(bolsig_pointer(66))
  rrt(074) = bolsig_rates(bolsig_pointer(67))
  rrt(075) = bolsig_rates(bolsig_pointer(68))
  rrt(076) = bolsig_rates(bolsig_pointer(69))
  rrt(077) = bolsig_rates(bolsig_pointer(70))
  rrt(078) = bolsig_rates(bolsig_pointer(71))
  rrt(079) = bolsig_rates(bolsig_pointer(72))
  rrt(080) = bolsig_rates(bolsig_pointer(73))
  rrt(081) = bolsig_rates(bolsig_pointer(74))
  rrt(082) = bolsig_rates(bolsig_pointer(75))
  rrt(083) = bolsig_rates(bolsig_pointer(76))
  rrt(084) = bolsig_rates(bolsig_pointer(77))
  rrt(085) = bolsig_rates(bolsig_pointer(78))
  rrt(086) = bolsig_rates(bolsig_pointer(79))
  rrt(087) = bolsig_rates(bolsig_pointer(80))
  rrt(088) = bolsig_rates(bolsig_pointer(81))
  rrt(089) = bolsig_rates(bolsig_pointer(82))
  rrt(090) = bolsig_rates(bolsig_pointer(83))
  rrt(091) = bolsig_rates(bolsig_pointer(84))
  rrt(092) = bolsig_rates(bolsig_pointer(85))
  rrt(093) = bolsig_rates(bolsig_pointer(86))
  rrt(094) = bolsig_rates(bolsig_pointer(87))
  rrt(095) = bolsig_rates(bolsig_pointer(88))
  rrt(096) = bolsig_rates(bolsig_pointer(89))
  rrt(097) = bolsig_rates(bolsig_pointer(90))
  rrt(098) = bolsig_rates(bolsig_pointer(91))
  rrt(099) = bolsig_rates(bolsig_pointer(92))
  rrt(100) = bolsig_rates(bolsig_pointer(93))
  rrt(101) = bolsig_rates(bolsig_pointer(94))
  rrt(102) = bolsig_rates(bolsig_pointer(95))
  rrt(103) = bolsig_rates(bolsig_pointer(96))
  rrt(104) = bolsig_rates(bolsig_pointer(97))
  rrt(105) = bolsig_rates(bolsig_pointer(98))
  rrt(106) = bolsig_rates(bolsig_pointer(99))
  rrt(107) = bolsig_rates(bolsig_pointer(100))
  rrt(108) = bolsig_rates(bolsig_pointer(101))
  rrt(109) = bolsig_rates(bolsig_pointer(102))
  rrt(110) = bolsig_rates(bolsig_pointer(103))
  rrt(111) = bolsig_rates(bolsig_pointer(104))
  rrt(112) = bolsig_rates(bolsig_pointer(105))
  rrt(113) = bolsig_rates(bolsig_pointer(106))
  rrt(114) = bolsig_rates(bolsig_pointer(107))
  rrt(115) = bolsig_rates(bolsig_pointer(108))
  rrt(116) = bolsig_rates(bolsig_pointer(109))
  rrt(117) = bolsig_rates(bolsig_pointer(110))
  rrt(118) = bolsig_rates(bolsig_pointer(111))
  rrt(119) = bolsig_rates(bolsig_pointer(112))
  rrt(120) = bolsig_rates(bolsig_pointer(113))
  rrt(121) = bolsig_rates(bolsig_pointer(114))
  rrt(122) = 2.6E-4
  rrt(123) = 1.1E1
  rrt(124) = bolsig_rates(bolsig_pointer(115))
  rrt(125) = bolsig_rates(bolsig_pointer(116))
  rrt(126) = bolsig_rates(bolsig_pointer(117))
  rrt(127) = bolsig_rates(bolsig_pointer(118))
  rrt(128) = bolsig_rates(bolsig_pointer(119))
  rrt(129) = bolsig_rates(bolsig_pointer(120))
  rrt(130) = bolsig_rates(bolsig_pointer(121))
  rrt(131) = bolsig_rates(bolsig_pointer(122))
  rrt(132) = bolsig_rates(bolsig_pointer(123))
  rrt(133) = bolsig_rates(bolsig_pointer(124))
  rrt(134) = bolsig_rates(bolsig_pointer(125))
  rrt(135) = bolsig_rates(bolsig_pointer(126))
  rrt(136) = bolsig_rates(bolsig_pointer(127))
  rrt(137) = bolsig_rates(bolsig_pointer(128))
  rrt(138) = bolsig_rates(bolsig_pointer(129))
  rrt(139) = bolsig_rates(bolsig_pointer(130))
  rrt(140) = bolsig_rates(bolsig_pointer(131))
  rrt(141) = bolsig_rates(bolsig_pointer(132))
  rrt(142) = bolsig_rates(bolsig_pointer(133))
  rrt(143) = bolsig_rates(bolsig_pointer(134))
  rrt(144) = bolsig_rates(bolsig_pointer(135))
  rrt(145) = bolsig_rates(bolsig_pointer(136))
  rrt(146) = bolsig_rates(bolsig_pointer(137))
  rrt(147) = bolsig_rates(bolsig_pointer(138))
  rrt(148) = bolsig_rates(bolsig_pointer(139))
  rrt(149) = bolsig_rates(bolsig_pointer(140))
  rrt(150) = bolsig_rates(bolsig_pointer(141))
  rrt(151) = bolsig_rates(bolsig_pointer(142))
  rrt(152) = bolsig_rates(bolsig_pointer(143))
  rrt(153) = bolsig_rates(bolsig_pointer(144))
  rrt(154) = bolsig_rates(bolsig_pointer(145))
  rrt(155) = bolsig_rates(bolsig_pointer(146))
  rrt(156) = bolsig_rates(bolsig_pointer(147))
  rrt(157) = 2.0E-05*TE**(-0.50)/TGAS
  rrt(158) = 3.939E-07*TE**(-0.40)
  rrt(159) = 1.608E-07*TE**(-0.50)
  rrt(160) = 3.683E-08*TE**(-0.55)
  rrt(161) = 9.64E-06*TE**(-0.34)
  rrt(162) = 5.4E-08*TE**(-0.70)
  rrt(163) = rrt(157)
  rrt(164) = 1.79E-08*TE**(-0.50)
  rrt(165) = 3E-30*ANY_NEUTRAL
  rrt(166) = 5E-31*TE**(-0.50)*ANY_NEUTRAL
  rrt(167) = 1.00E-31*ANY_NEUTRAL
  rrt(168) = 1.00E-26*ANY_NEUTRAL
  rrt(169) = 6E-07*TE**(-0.50)*TGAS**(-0.50)
  rrt(170) = rrt(168)
  rrt(171) = 2.251E-07*TE**(-0.5)
  rrt(172) = 2.8E-11*EXP(-26500/TGAS)
  rrt(173) = 1.00E-15
  rrt(174) = 2.80E-11*EXP(-MAX(0.0,26500-0.45*9.60E+02)/TGAS)
  rrt(175) = 2.80E-11*EXP(-MAX(0.0,26500-0.45*1.92E+03)/TGAS)
  rrt(176) = 2.80E-11*EXP(-MAX(0.0,26500-0.45*3.38E+03)/TGAS)
  rrt(177) = 2.80E-11*EXP(-MAX(0.0,26500-0.45*3.85E+03)/TGAS)
  rrt(178) = 1E-15
  rrt(179) = 1E-15
  rrt(180) = 1E-15
  rrt(181) = 1E-15
  rrt(182) = 8.2E-34*EXP(-1510/TGAS)*ANY_NEUTRAL
  rrt(183) = 4.00E-25
  rrt(184) = 6.5E-32*ANY_NEUTRAL
  rrt(185) = 4.2E-12*EXP(-24000/TGAS)
  rrt(186) = 8.20E-34*EXP(-1510.0/TGAS)*ANY_NEUTRAL
  rrt(187) = 4.20E-12*EXP(-MAX(0.0,24000-0.45*2.24E+03)/TGAS)
  rrt(188) = 4.20E-12*EXP(-MAX(0.0,24000-0.45*4.44E+03)/TGAS)
  rrt(189) = 4.20E-12*EXP(-MAX(0.0,24000-0.45*6.62E+03)/TGAS)
  rrt(190) = 4.20E-12*EXP(-MAX(0.0,24000-0.45*0.00E+00)/TGAS)
  rrt(191) = rrt(187)
  rrt(192) = rrt(188)
  rrt(193) = rrt(189)
  rrt(194) = 4.0E-25
  rrt(195) = 6.50E-32*ANY_NEUTRAL
  rrt(196) = 3.0E-11
  rrt(197) = 3.3E-13
  rrt(198) = 1.81E-33*(TGAS/300)**(-1.20)*ANY_NEUTRAL
  rrt(199) = 3.00E-11
  rrt(200) = 3.00E-11
  rrt(201) = 3.00E-11
  rrt(202) = 1.81E-33*((TGAS/300.0)**(-1.20))*ANY_NEUTRAL
  rrt(203) = rrt(202)
  rrt(204) = rrt(202)
  rrt(205) = 2.136E-29*(TGAS/300)**(-3.08)*EXP(-2114/TGAS)*ANY_NEUTRAL
  rrt(206) = 5.0E-11
  rrt(207) = 3.1E-14*TGAS**(0.75)*EXP(-1575/TGAS)
  rrt(208) = 4.1175E-10*EXP(-11430/TGAS)*ANY_NEUTRAL
  rrt(209) = 1.27E-32*EXP(-170/TGAS)/(TGAS/300)*ANY_NEUTRAL
  rrt(210) = 2.3E-29*ANY_NEUTRAL
  rrt(211) = 9.4E-10
  rrt(212) = 4.5E-10
  rrt(213) = 1.1E-09
  rrt(214) = 1.0E-09
  rrt(215) = 9.0E-29*ANY_NEUTRAL
  rrt(216) = 1.0E-29*ANY_NEUTRAL
  rrt(217) = 5.5E-10
  rrt(218) = 4.8E-10
  rrt(219) = 3.0E-28*ANY_NEUTRAL
  rrt(220) = 9.40E-10
  rrt(221) = 9.40E-10
  rrt(222) = 9.40E-10
  rrt(223) = 9.40E-10
  rrt(224) = 4.50E-10
  rrt(225) = 4.50E-10
  rrt(226) = 4.50E-10
  rrt(227) = 4.50E-10
  rrt(228) = 1.10E-09
  rrt(229) = 1.10E-09
  rrt(230) = 1.10E-09
  rrt(231) = 1.10E-09
  rrt(232) = 1.00E-09
  rrt(233) = 1.00E-09
  rrt(234) = 1.00E-09
  rrt(235) = 1.00E-09
  rrt(236) = 9.00E-29*ANY_NEUTRAL
  rrt(237) = rrt(236)
  rrt(238) = rrt(236)
  rrt(239) = rrt(236)
  rrt(240) = 1.00E-29*ANY_NEUTRAL
  rrt(241) = rrt(240)
  rrt(242) = rrt(240)
  rrt(243) = rrt(240)
  rrt(244) = 5.50E-10
  rrt(245) = 5.50E-10
  rrt(246) = 5.50E-10
  rrt(247) = 5.50E-10
  rrt(248) = 4.80E-10
  rrt(249) = 4.80E-10
  rrt(250) = 4.80E-10
  rrt(251) = 4.80E-10
  rrt(252) = 3.00E-28*ANY_NEUTRAL
  rrt(253) = rrt(252)
  rrt(254) = rrt(252)
  rrt(255) = rrt(252)
  rrt(256) = 2.30E-29*ANY_NEUTRAL
  rrt(257) = rrt(256)
  rrt(258) = rrt(256)
  rrt(259) = rrt(256)
  rrt(260) = 9.40E-10*(12.06/6.23)**2
  rrt(261) = 9.40E-10*(12.06/11.1)**2
  rrt(262) = 4.50E-10*(13.8/6.23)**2
  rrt(263) = 4.50E-10*(13.8/11.1)**2
  rrt(264) = 1.10E-09*(17.51/6.23)**2
  rrt(265) = 1.10E-09*(17.51/11.1)**2
  rrt(266) = 1.00E-09*(13.8/6.23)**2
  rrt(267) = 1.00E-09*(13.8/11.1)**2
  rrt(268) = rrt(236)
  rrt(269) = rrt(236)
  rrt(270) = rrt(240)
  rrt(271) = rrt(240)
  rrt(272) = 5.50E-10
  rrt(273) = 5.50E-10
  rrt(274) = 4.80E-10
  rrt(275) = 4.80E-10
  rrt(276) = rrt(252)
  rrt(277) = rrt(252)
  rrt(278) = rrt(256)
  rrt(279) = rrt(256)
  rrt(280) = rrt(256)
  rrt(281) = rrt(256)
  rrt(282) = 4.9E-12*(TGAS/300)**(0.50)*EXP(-4580/TGAS)
  rrt(283) = 5.5E-10
  rrt(284) = 5.0E-13
  rrt(285) = 1.1E-09
  rrt(286) = 9.0E-10
  rrt(287) = 2.6E-26*ANY_NEUTRAL
  rrt(288) = 4.2E-26*ANY_NEUTRAL
  rrt(289) = 5.0E-13
  rrt(290) = 4.90E-12*((TGAS/300.0)**(0.5))*EXP(-4580.0/TGAS)
  rrt(291) = 5.50E-10
  rrt(292) = 5.0E-13
  rrt(293) = 1.10E-09
  rrt(294) = 9.00E-10
  rrt(295) = 2.60E-26*ANY_NEUTRAL
  rrt(296) = 4.20E-26*ANY_NEUTRAL
  rrt(297) = 5.00E-13
  rrt(298) = 5.50E-10
  rrt(299) = 5.50E-10
  rrt(300) = 5.50E-10
  rrt(301) = 5.50E-10
  rrt(302) = 5.00E-13
  rrt(303) = 5.00E-13
  rrt(304) = 5.00E-13
  rrt(305) = 5.00E-13
  rrt(306) = 1.10E-09
  rrt(307) = 1.10E-09
  rrt(308) = 1.10E-09
  rrt(309) = 1.10E-09
  rrt(310) = 9.00E-10
  rrt(311) = 9.00E-10
  rrt(312) = 9.00E-10
  rrt(313) = 9.00E-10
  rrt(314) = rrt(295)
  rrt(315) = rrt(295)
  rrt(316) = rrt(295)
  rrt(317) = rrt(295)
  rrt(318) = rrt(296)
  rrt(319) = rrt(296)
  rrt(320) = rrt(296)
  rrt(321) = rrt(296)
  rrt(322) = 5.00E-13*(17.51/6.22)**2
  rrt(323) = 5.00E-13*(17.51/7.9)**2
  rrt(324) = 5.00E-13*(17.51/13.5)**2
  rrt(325) = 5.00E-13*(17.51/10.01)**2
  rrt(326) = 5.3E-11
  rrt(327) = 1.2E-10
  rrt(328) = 5.0E-12
  rrt(329) = 6.2E-10
  rrt(330) = 3.8E-10
  rrt(331) = 1.9E-10*(TGAS/300)**(-0.50)
  rrt(332) = 2.4E-30*ANY_NEUTRAL
  rrt(333) = 3.5E-31*ANY_NEUTRAL
  rrt(334) = 1.00E-12
  rrt(335) = 3.0E-28/(TGAS/300)*ANY_NEUTRAL
  rrt(336) = 2.18E-18
  rrt(337) = 2.3E-11
  rrt(338) = 5.30E-11
  rrt(339) = 5.30E-11
  rrt(340) = 5.30E-11
  rrt(341) = 1.20E-10
  rrt(342) = 1.20E-10
  rrt(343) = 1.20E-10
  rrt(344) = 5.00E-12
  rrt(345) = 5.00E-12
  rrt(346) = 5.00E-12
  rrt(347) = 6.20E-10
  rrt(348) = 6.20E-10
  rrt(349) = 6.20E-10
  rrt(350) = 3.80E-10
  rrt(351) = 3.80E-10
  rrt(352) = 3.80E-10
  rrt(353) = 1.90E-11*((TGAS/300.0)**(-0.50))
  rrt(354) = rrt(353)
  rrt(355) = rrt(353)
  rrt(356) = 2.40E-30*ANY_NEUTRAL
  rrt(357) = rrt(356)
  rrt(358) = rrt(356)
  rrt(359) = 3.50E-31*ANY_NEUTRAL
  rrt(360) = rrt(359)
  rrt(361) = rrt(359)
  rrt(362) = 1.00E-12
  rrt(363) = 1.00E-12
  rrt(364) = 1.00E-12
  rrt(365) = 3.00E-28*((TGAS/300.0)**(-1.0))*ANY_NEUTRAL
  rrt(366) = rrt(365)
  rrt(367) = rrt(365)
  rrt(368) = 2.18E-18
  rrt(369) = 2.18E-18
  rrt(370) = 2.18E-18
  rrt(371) = 2.30E-11
  rrt(372) = 2.30E-11
  rrt(373) = 2.30E-11
  rrt(374) = 5.30E-11*(12.06/0.98)**2
  rrt(375) = 5.30E-11*(12.06/8.4)**2
  rrt(376) = 1.20E-10*(12.06/0.98)**2
  rrt(377) = 1.20E-10*(12.06/8.4)**2
  rrt(378) = 5.00E-12
  rrt(379) = 5.00E-12
  rrt(380) = 6.20E-10
  rrt(381) = 6.20E-10
  rrt(382) = 3.80E-10
  rrt(383) = 3.80E-10
  rrt(384) = 1.90E-11*((TGAS/300.0)**(-0.50))*(12.06/0.98)**2
  rrt(385) = 1.90E-11*((TGAS/300.0)**(-0.50))*(12.06/8.4)**2
  rrt(386) = rrt(356)
  rrt(387) = rrt(356)
  rrt(388) = rrt(359)
  rrt(389) = rrt(359)
  rrt(390) = 3.00E-10
  rrt(391) = 3.00E-10
  rrt(392) = rrt(365)
  rrt(393) = rrt(365)
  rrt(394) = 2.18E-18
  rrt(395) = 2.18E-18
  rrt(396) = 2.30E-11
  rrt(397) = 2.30E-11
  rrt(398) = 1.1E-10
  rrt(399) = 5.2E-11
  rrt(400) = 5.2E-11
  rrt(401) = 1.1E-10
  rrt(402) = 1.64E-10
  rrt(403) = 9.62E-11
  rrt(404) = 5.0E-07
  rrt(405) = 5.0E-07
  rrt(406) = 6.0E-07
  rrt(407) = 1.4E-10
  rrt(408) = 1.0E-12*ANY_NEUTRAL
  rrt(409) = 5.0E-07
  rrt(410) = 5.0E-07
  rrt(411) = 6.0E-07
  rrt(412) = 5.0E-07
  rrt(413) = 5.0E-07
  rrt(414) = 6.0E-07
  rrt(415) = 1.0E-14*ANY_NEUTRAL
  rrt(416) = 5.0E-07
  rrt(417) = 5.0E-07
  rrt(418) = 6.0E-07
  rrt(419) = 3.0E-07
  rrt(420) = 8.0E-11
  rrt(421) = 3.0E-07
  rrt(422) = 1.1E-10
  rrt(423) = 1.4E-11
  rrt(424) = 1.4E-11
  rrt(425) = 1.3E-10
  rrt(426) = 8.0E-10
  rrt(427) = 3.0E-10
  rrt(428) = 4.0E-10
  rrt(429) = 3.0E-10
  rrt(430) = 1.00E-10
  rrt(431) = rrt(240)
  rrt(432) = 2.3E-10
  rrt(433) = 3.3E-10
  rrt(434) = 3.3E-10
  rrt(435) = 1.00E-13
  rrt(436) = 1.00E-13
  rrt(437) = 2.5E-10
  rrt(438) = 4.0E-10
  rrt(439) = 3.0E-10
  rrt(440) = 3.0E-10
  rrt(441) = 2.0E-25*ANY_NEUTRAL
  rrt(442) = 2.7E-07
  rrt(443) = 2.0E-07
  rrt(444) = 4.2E-07
  rrt(445) = rrt(441)
  rrt(446) = 2.7E-10*(TGAS/300)**(0.50)*EXP(-5590/TGAS)*ANY_NEUTRAL
  rrt(447) = 2.0E-07
  rrt(448) = 1.0E-07
  rrt(449) = 1.0E-07
  rrt(450) = 2.3E-11*ANY_NEUTRAL
  rrt(451) = 4.0E-08
  rrt(452) = rrt(441)
  rrt(453) = 1.0E-07
  rrt(454) = 2.6E-08
  rrt(455) = rrt(441)
  rrt(456) = 4.0E-12*ANY_NEUTRAL
  rrt(457) = 3.08E-12*ANY_NEUTRAL
  rrt(458) = 1.73E-13*ANY_NEUTRAL
  rrt(459) = 1.07E-14
  rrt(460) = 7.48E-15
  rrt(461) = 7.48E-15
  rrt(462) = 9.00E-18
  rrt(463) = 2.79E-17
  rrt(464) = 2.79E-17
  rrt(465) = 2.90E-14
  rrt(466) = 2.03E-14
  rrt(467) = 2.03E-14
  rrt(468) = 7.72E-16
  rrt(469) = 2.32E-16
  rrt(470) = 3.09E-16
  rrt(471) = 6.05E-15
  rrt(472) = 1.81E-15
  rrt(473) = 2.42E-15
  rrt(474) = 2.42E-15
  rrt(475) = 1.70E-18
  rrt(476) = 5.10E-19
  rrt(477) = 6.80E-19
  rrt(478) = 4.33E-14
  rrt(479) = 3.03E-14
  rrt(480) = 3.03E-14
  rrt(481) = 9.08E-18
  rrt(482) = 6.18E-15
  rrt(483) = 6.18E-15
  rrt(484) = 1.34E-23
  rrt(485) = 1.34E-23
  rrt(486) = 4.78E-24
  rrt(487) = 7.55E-23
  rrt(488) = 2.52E-23
  rrt(489) = 2.52E-23
  rrt(490) = 7.55E-23
  rrt(491) = 2.52E-23
  rrt(492) = 2.52E-23
  rrt(493) = 7.55E-23
  rrt(494) = 2.52E-23
  rrt(495) = 2.52E-23
  rrt(496) = 1.0E-11*ANY_NEUTRAL
  rrt(497) = rrt(496)
  rrt(498) = 1.2E-11*ANY_NEUTRAL
  rrt(499) = rrt(498)
  rrt(500) = rrt(498)
  rrt(501) = rrt(498)
  rrt(502) = 3.8E-18*EXP(-205.0/TGAS)*ANY_NEUTRAL
  rrt(503) = 3.0E-13*ANY_NEUTRAL
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
