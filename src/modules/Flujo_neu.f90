!================================================================================
! MÓDULO: Flujo_neu_M2
!   Espectro combinado Kopeikin (E<2 MeV) + Mueller (E>=2 MeV)
!   Normalizado al flujo total integrado total_flux_cm2_s = 1.5e13 cm⁻² s⁻¹
!================================================================================
module Flujo_neu
  use constan
  implicit none

  ! ---------- Datos de Kopeikin (50 puntos) ----------
  integer, parameter :: n_kop = 50
  real(dp), parameter :: E_kop(n_kop) = (/ &
       0.010_dp, 0.020_dp, 0.035_dp, 0.040_dp, 0.070_dp, 0.100_dp, 0.130_dp, 0.160_dp, &
       0.165_dp, 0.180_dp, 0.215_dp, 0.230_dp, 0.280_dp, 0.330_dp, 0.335_dp, 0.350_dp, &
       0.390_dp, 0.400_dp, 0.435_dp, 0.440_dp, 0.500_dp, 0.700_dp, 0.900_dp, 1.000_dp, &
       1.185_dp, 1.190_dp, 1.250_dp, 1.300_dp, 1.500_dp, 1.700_dp, 1.800_dp, 1.900_dp, &
       2.000_dp, 2.250_dp, 2.500_dp, 2.750_dp, 3.000_dp, 3.250_dp, 3.500_dp, 4.000_dp, &
       4.500_dp, 5.000_dp, 5.500_dp, 6.000_dp, 6.500_dp, 7.000_dp, 7.500_dp, 8.000_dp, &
       8.500_dp, 9.000_dp /)
  
  real(dp), parameter :: rho_kop(n_kop) = (/ &
       0.774e-1_dp, 0.301_dp, 0.848_dp, 0.354_dp, 0.989_dp, 1.81_dp, 2.74_dp, 3.85_dp, &
       3.26_dp, 3.74_dp, 4.66_dp, 4.08_dp, 5.03_dp, 5.97_dp, 4.30_dp, 4.08_dp, &
       4.42_dp, 4.04_dp, 4.39_dp, 2.84_dp, 3.03_dp, 3.14_dp, 3.23_dp, 3.11_dp, &
       2.80_dp, 2.24_dp, 2.06_dp, 1.77_dp, 1.59_dp, 1.48_dp, 1.42_dp, 1.36_dp, &
       1.30_dp, 1.08_dp, 0.882_dp, 0.733_dp, 0.611_dp, 0.505_dp, 0.411_dp, 0.261_dp, &
       0.156_dp, 0.928e-1_dp, 0.549e-1_dp, 0.317e-1_dp, 0.174e-1_dp, 0.897e-2_dp, &
       0.366e-2_dp, 0.119e-2_dp, 0.299e-3_dp, 0.831e-4_dp /)

  ! ---------- Coeficientes de Mueller ----------
  integer, parameter :: n_iso = 4
  real(dp), parameter :: fission_frac(n_iso) = (/ 0.58_dp, 0.07_dp, 0.30_dp, 0.05_dp /)
  
  real(dp), parameter :: U235(6) = (/ 3.217_dp, -3.111_dp, 1.395_dp, -0.369_dp, 0.04445_dp, -0.002053_dp /)
  real(dp), parameter :: U238(6) = (/ 0.4833_dp, 0.1927_dp, -0.1283_dp, -0.006762_dp, 0.002233_dp, -0.0001536_dp /)
  real(dp), parameter :: Pu239(6)= (/ 6.413_dp, -7.432_dp, 3.535_dp, -0.882_dp, 0.1025_dp, -0.00455_dp /)
  real(dp), parameter :: Pu241(6)= (/ 3.251_dp, -3.204_dp, 1.428_dp, -0.3675_dp, 0.04254_dp, -0.001896_dp /)

  real(dp) :: spectrum_integral = 0.0_dp
  logical :: norm_initialized = .false.

contains

  ! Función pura para asegurar que no hay efectos secundarios durante la interpolación
  pure function kopeikin_spectrum(E_nu) result(rho)
    real(dp), intent(in) :: E_nu
    real(dp) :: rho
    integer :: i
    
    if (E_nu < E_kop(1) .or. E_nu > E_kop(n_kop)) then
       rho = 0.0_dp
       return
    end if
    
    do i = 1, n_kop-1
       if (E_nu >= E_kop(i) .and. E_nu <= E_kop(i+1)) then
          rho = rho_kop(i) + (rho_kop(i+1) - rho_kop(i)) * (E_nu - E_kop(i)) / (E_kop(i+1) - E_kop(i))
          return
       end if
    end do
    rho = rho_kop(n_kop)
  end function kopeikin_spectrum

  pure function mueller_spectrum(E_nu, iso) result(s)
    real(dp), intent(in) :: E_nu
    integer, intent(in) :: iso
    real(dp) :: s, log_s
    
    select case(iso)
    case(1)
       log_s = U235(1) + U235(2)*E_nu + U235(3)*(E_nu**2) + U235(4)*(E_nu**3) &
             + U235(5)*(E_nu**4) + U235(6)*(E_nu**5)
    case(2)
       log_s = U238(1) + U238(2)*E_nu + U238(3)*(E_nu**2) + U238(4)*(E_nu**3) &
             + U238(5)*(E_nu**4) + U238(6)*(E_nu**5)
    case(3)
       log_s = Pu239(1) + Pu239(2)*E_nu + Pu239(3)*(E_nu**2) + Pu239(4)*(E_nu**3) &
             + Pu239(5)*(E_nu**4) + Pu239(6)*(E_nu**5)
    case(4)
       log_s = Pu241(1) + Pu241(2)*E_nu + Pu241(3)*(E_nu**2) + Pu241(4)*(E_nu**3) &
             + Pu241(5)*(E_nu**4) + Pu241(6)*(E_nu**5)
    end select
    s = exp(log_s)
  end function mueller_spectrum

  pure function espectro_total(E_nu) result(rho)
    real(dp), intent(in) :: E_nu
    real(dp) :: rho
    integer :: i
    
    if (E_nu < 0.0_dp .or. E_nu > Enu_max) then
       rho = 0.0_dp
       return
    end if
    
    if (E_nu < 2.0_dp) then
       rho = kopeikin_spectrum(E_nu)
    else
       rho = 0.0_dp
       do i = 1, n_iso
          rho = rho + fission_frac(i) * mueller_spectrum(E_nu, i)
       end do
    end if
  end function espectro_total

  ! Integración numérica optimizada con casteo explícito
  function compute_spectrum_integral() result(intg)
    real(dp) :: intg
    integer, parameter :: n_int = 10000 ! Aumentado para mayor precisión inicial
    real(dp) :: dE, E
    integer :: i
    
    dE = Enu_max / real(n_int - 1, dp)
    intg = 0.5_dp * (espectro_total(0.0_dp) + espectro_total(Enu_max))
    
    do i = 2, n_int-1
       E = real(i-1, dp) * dE
       intg = intg + espectro_total(E)
    end do
    intg = intg * dE
  end function compute_spectrum_integral

  subroutine Initialize_Flux()
    if (.not. norm_initialized) then
       spectrum_integral = compute_spectrum_integral()
       norm_initialized = .true.
       print *, "=== ESPECTRO COMBINADO (Kopeikin + Mueller) ==="
       print *, "Integral del espectro (0-10 MeV):", spectrum_integral
       print *, "Factor de normalización:", total_flux_cm2_s / spectrum_integral
       print *, "==============================================="
    end if
  end subroutine Initialize_Flux

  ! Función para invocar el flujo final normalizado
  function Flux(Enu) result(flujo)
    real(dp), intent(in) :: Enu
    real(dp) :: flujo
    
    if (.not. norm_initialized) then
       print *, "ERROR: Initialize_Flux no ha sido llamado antes de Flux."
       stop
    end if
    flujo = total_flux_cm2_s * espectro_total(Enu) / spectrum_integral
  end function Flux

end module Flujo_neu
