module constan
  implicit none
  ! Definición de precisión estandarizada (mínimo 15 decimales, exponente hasta 307)
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! Constantes fundamentales
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: GF_MeV = 1.1663787e-11_dp
  real(dp), parameter :: alpha_EM = 1.0_dp / 137.035999084_dp
  real(dp), parameter :: sin2W = 0.23857_dp
  
  ! Masas (MeV)
  real(dp), parameter :: m_e_MeV = 0.5109989461_dp
  real(dp), parameter :: m_n_MeV = 939.565420_dp
  real(dp), parameter :: m_p_MeV = 938.272088_dp
  real(dp), parameter :: u_MeV = 931.49410242_dp

  !=====================================================================
  ! PARÁMETROS BSM: Fermiones Estériles y Mediador Escalar
  ! Masas de los quarks ligeros (MeV) - Valores estándar PDG
  real(dp), parameter :: m_u_MeV = 2.16_dp
  real(dp), parameter :: m_d_MeV = 4.67_dp
  
  ! Fracciones de masa del quark en el nucleón (Paper Ec. 11, pag 5)
  real(dp), parameter :: f_Tu_p = 0.026_dp
  real(dp), parameter :: f_Td_p = 0.038_dp
  real(dp), parameter :: f_Tu_n = 0.018_dp
  real(dp), parameter :: f_Td_n = 0.056_dp
  !=====================================================================
  
  ! Factores de conversión
  real(dp), parameter :: hbarc = 197.3269804e-13_dp
  real(dp), parameter :: conv_MeV2_to_cm2 = (hbarc)**2   ! 3.8938e-22 cm² MeV²
  real(dp), parameter :: hbarc_fm = 197.3269804_dp

  ! Propiedades del blanco de Germanio (Promedio efectivo para CEνNS)
  real(dp), parameter :: A_Ge = 72.6_dp
  real(dp), parameter :: Z_Ge = 32.0_dp
  real(dp), parameter :: N_Ge = A_Ge - Z_Ge
  real(dp), parameter :: M_nucleus_Ge = A_Ge * u_MeV

  ! Parámetros del detector CONUS+
  real(dp), parameter :: eta_Ge = 2.96e-6_dp          ! MeV 
  real(dp), parameter :: sigma0_eV = 20.38_dp
  real(dp), parameter :: sigma0 = sigma0_eV * 1.0e-6_dp   ! MeV
  real(dp), parameter :: Fano_factor = 0.1096_dp

  ! Parámetros de Quenching (Teoría de Lindhard)
  real(dp), parameter :: k_Lindhard = 0.162_dp
  real(dp), parameter :: Z_Lindhard = 32.0_dp

  ! Límites de integración y rangos (Enu_min_global se calcula dinámicamente)
  real(dp), parameter :: Enu_max = 15.0_dp            ! MeV
  real(dp), parameter :: Eer_min_global = 2.96e-6_dp  ! MeV
  real(dp), parameter :: Eer_max_global = 0.001_dp    ! MeV 

  ! Parámetros computacionales y experimentales
  integer, parameter :: nEnu_points = 100
  integer, parameter :: nEer_points = 100
  real(dp), parameter :: total_flux_cm2_s = 1.5e13_dp ! Flujo integrado total (cm⁻² s⁻¹)

  ! Constantes molares y de avogadro
  real(dp), parameter :: Avogadro = 6.02214076e23_dp
  real(dp), parameter :: M_Ge_gmol = 72.63_dp
  real(dp), parameter :: N_nuclei_per_kg = Avogadro / (M_Ge_gmol / 1000.0_dp) ! ~8.29e24
end module constan
