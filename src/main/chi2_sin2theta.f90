!================================================================================
! PROGRAMA: chi2_angle
! Calcula el ajuste de chi-cuadrado para determinar el ángulo de mezcla débil 
! (sin²θ_W) usando el espectro del detector efectivo de CONUS+ y minimizando 
! analíticamente el parámetro de ruido α.
!================================================================================
program chi2_sin2theta
  use constan
  use cross_sec
  use Flujo_neu
  use Fun_resolu
  use Integrales
  implicit none

  integer, parameter :: nbins = 19
  real(dp) :: bin_edges(nbins+1)
  real(dp) :: R_exp(nbins), sigma_exp(nbins)
  real(dp) :: R_ref(nbins), R_th(nbins)
  real(dp) :: s2w_min, s2w_max, s2w, qw, qw_ref, factor, chi2, alpha_best, sum1, sum2
  real(dp) :: sin2W_ref
  integer :: step, n_steps, i
  real(dp), parameter :: sigma_alpha = 0.169_dp
  character(100) :: filename

  print *, "===== PROGRAMA CHI2_ANGLE (para sin²θ_W) ====="

  ! Datos experimentales (excesos de eventos C3, digitalizados Fig 1)
  R_exp = [13.27102804_dp, 38.06853583_dp, 29.09657321_dp, 21.37071651_dp, 12.2741433_dp, &
           4.92211838_dp, 4.29906542_dp, 0.0623053_dp, 9.1588785_dp, 2.42990654_dp, &
           3.92523364_dp, -2.55451713_dp, 0.56074766_dp, -10.52959502_dp, -4.29906542_dp, &
           -8.78504673_dp, -0.56074766_dp, -0.93457944_dp, 13.52024922_dp]

  sigma_exp = [22.9906542_dp, 15.70093458_dp, 12.64485981_dp, 12.21183801_dp, 11.52647975_dp, &
               11.58878505_dp, 11.27725857_dp, 10.7788162_dp, 11.02803738_dp, 10.28037383_dp, &
               10.34267913_dp, 10.7788162_dp, 10.34267913_dp, 10.15576324_dp, 10.09345794_dp, &
               9.90654206_dp, 9.84423676_dp, 9.78193146_dp, 9.59501558_dp]

  call Define_Energy_Bins(nbins, bin_edges)
  call Initialize_Flux()

  ! Generar Espectro de Referencia con el valor SM
  sin2W_ref = 0.23857_dp
  call compute_R_ref(bin_edges, R_ref)

  ! Carga débil teórica de referencia
  qw_ref = (0.5_dp - 2.0_dp*sin2W_ref) * Z_Ge + (-0.5_dp) * N_Ge

  ! Barrido de parámetros
  s2w_min = 0.01_dp
  s2w_max = 0.50_dp
  n_steps = 3000
  filename = "chi2_sin2theta.dat"
  open(10, file=filename, status='replace')
  write(10, '(A)') "# sin2theta_W   chi2"

  print *, "Minimizando chi-cuadrado sobre", n_steps, "puntos..."

  do step = 0, n_steps
     ! Casteo explícito a dp en el bucle numérico
     s2w = s2w_min + real(step, dp) * (s2w_max - s2w_min) / real(n_steps, dp)
     qw = (0.5_dp - 2.0_dp*s2w) * Z_Ge + (-0.5_dp) * N_Ge
     factor = (qw / qw_ref)**2
     R_th = R_ref * factor

     ! Minimización analítica de alpha (Nuisance parameter)
     sum1 = 0.0_dp; sum2 = 0.0_dp
     do i = 1, nbins
        sum1 = sum1 + (R_exp(i) * R_th(i)) / (sigma_exp(i)**2)
        sum2 = sum2 + (R_th(i)**2) / (sigma_exp(i)**2)
     end do
     alpha_best = (sum1 - sum2) / (sum2 + 1.0_dp / sigma_alpha**2)

     ! Cálculo final de chi2
     chi2 = 0.0_dp
     do i = 1, nbins
        chi2 = chi2 + ((R_exp(i) - (1.0_dp + alpha_best) * R_th(i))**2) / (sigma_exp(i)**2)
     end do
     chi2 = chi2 + (alpha_best / sigma_alpha)**2

     write(10, '(F10.6, 2X, F12.4)') s2w, chi2
  end do

  close(10)
  print *, "=== ANÁLISIS FINALIZADO ==="
  print *, "Archivo de datos generado: ", trim(filename)

contains

  subroutine compute_R_ref(bin_edges, R_ref_out)
    real(dp), intent(in) :: bin_edges(:)
    real(dp), intent(out) :: R_ref_out(:)
    integer :: iBin, j, k
    real(dp) :: E_reco_min, E_reco_max
    real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
    real(dp) :: sum_Eer, sum_Enu, gauss_int, dsigma, flux_val
    real(dp) :: Eer, E_nu, E_nu_min
    real(dp) :: exposure_kg_s
    
    ! Propiedades del Detector Efectivo (Figura 1)
    real(dp), parameter :: detector_mass_kg = 1.0_dp
    real(dp), parameter :: detector_days = 119.0_dp

    exposure_kg_s = detector_mass_kg * detector_days * 86400.0_dp
    R_ref_out = 0.0_dp

    do iBin = 1, size(bin_edges)-1
       E_reco_min = bin_edges(iBin)
       E_reco_max = bin_edges(iBin+1)
       
       ! Integración de energía verdadera respetando el umbral físico del Ge
       Eer_min_int = Eer_min_global    
       Eer_max_int = E_reco_max + 5.0_dp * sigma0
       if (Eer_min_int >= Eer_max_int) cycle

       dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
       sum_Eer = 0.0_dp

       do j = 1, nEer_points
          Eer = Eer_min_int + real(j-1, dp) * dEer
          
          gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
          if (gauss_int < 1.0e-10_dp) cycle

          E_nu_min = sqrt(M_nucleus_Ge * T_from_Eer_NR(Eer)) / 2.0_dp
          if (E_nu_min > Enu_max) cycle

          dEnu = (Enu_max - E_nu_min) / real(nEnu_points - 1, dp)
          sum_Enu = 0.0_dp

          do k = 1, nEnu_points
             E_nu = E_nu_min + real(k-1, dp) * dEnu
             if (E_nu > Enu_max) exit

             dsigma = dSigmadEer_SM(Eer, E_nu)
             flux_val = Flux(E_nu)

             if (k == 1 .or. k == nEnu_points) then
                sum_Enu = sum_Enu + 0.5_dp * dsigma * flux_val
             else
                sum_Enu = sum_Enu + dsigma * flux_val
             end if
          end do
          sum_Enu = sum_Enu * dEnu * N_nuclei_per_kg

          if (j == 1 .or. j == nEer_points) then
             sum_Eer = sum_Eer + 0.5_dp * sum_Enu * gauss_int * dEer
          else
             sum_Eer = sum_Eer + sum_Enu * gauss_int * dEer
          end if
       end do

       R_ref_out(iBin) = sum_Eer * exposure_kg_s
    end do
  end subroutine compute_R_ref

end program chi2_sin2theta
