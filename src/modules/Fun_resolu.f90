!================================================================================
! MÓDULO: Fun_resolu_M
! Contiene: Función de resolución de los detectores CONUS+ y su integral analítica.
! Actualizado al estándar de precisión 'dp' y a la Ec. 30 del artículo.
!================================================================================
module Fun_resolu
  use constan
  implicit none
  
contains

  !----------------------------------------------------------------------------
  ! Cálculo del ancho de resolución (sigma_res) en función de Eer (Ecuación 30)
  ! Devuelve sigma en MeV
  !----------------------------------------------------------------------------
  pure real(dp) function Sigma_res(Eer)
    real(dp), intent(in) :: Eer
    
    ! Si la energía es negativa (físicamente imposible, pero seguro numérico)
    if (Eer <= 0.0_dp) then
       Sigma_res = sigma0 ! Solo ruido electrónico
       return
    end if
    
    ! Ec 30: sigma = sqrt(sigma0^2 + Fano * eta * Eer)
    Sigma_res = sqrt(sigma0**2 + Fano_factor * eta_Ge * Eer)
  end function Sigma_res

  !----------------------------------------------------------------------------
  ! Función de resolución gaussiana (Ecuación 29)
  ! Recibe energía reconstruida y energía verdadera (ionización)
  !----------------------------------------------------------------------------
  pure real(dp) function Gaussian(E_reco, Eer)
    real(dp), intent(in) :: E_reco, Eer
    real(dp) :: sigma, arg
    
    sigma = Sigma_res(Eer)
    arg = (E_reco - Eer) / sigma
    
    Gaussian = exp(-0.5_dp * arg**2) / (sigma * sqrt(2.0_dp * pi))
  end function Gaussian
  
!----------------------------------------------------------------------------
  ! Integral analítica de la gaussiana en el bin de energía [a, b]
  ! Devuelve la fracción de la señal en Eer que cae dentro del bin
  !----------------------------------------------------------------------------
  pure real(dp) function Integral_Gaussiana(a, b, Eer)
    real(dp), intent(in) :: a, b, Eer
    real(dp) :: sigma, za, zb, erf_za, erf_zb
    real(dp) :: sqrt2_inv
    
    sigma = Sigma_res(Eer)
    
    ! Constante correcta para multiplicar el argumento
    sqrt2_inv = 1.0_dp / sqrt(2.0_dp)
    
    ! Limites de integración normalizados ¡CORREGIDOS!
    ! Se multiplica por sqrt2_inv, lo que equivale a dividir por sqrt(2)
    za = (a - Eer) * sqrt2_inv / sigma
    zb = (b - Eer) * sqrt2_inv / sigma
    
    ! Función error intrínseca
    erf_za = erf(za)
    erf_zb = erf(zb)
    
    Integral_Gaussiana = 0.5_dp * (erf_zb - erf_za)
  end function Integral_Gaussiana
end module Fun_resolu
