!================================================================================
! PROGRAMA: chi2_mchi_scan.f90
! Objetivo: Escanear la masa del fermión estéril (m_chi) fijando el mediador (M_S)
!           y calcular el valor de Chi^2 con la incertidumbre sistemática (alpha).
!================================================================================
program chi2_mchi_scan
  use constan
  use cross_sec
  use Flujo_neu     
  use Fun_resolu
  use Integrales
  implicit none
   
  integer, parameter :: nBins = 19
  real(dp) :: bin_edges(nBins+1)
  real(dp) :: E_reco_min, E_reco_max, E_reco_center
  integer  :: iBin, j, k, step
  real(dp) :: Eer, E_nu, T_N, Enu_min_nuc, Enu_min_e, E_nu_min_loop
  real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
  real(dp) :: sum_Enu, sum_Eer, gauss_int, flux_val
  
  ! Arreglos para el análisis estadístico
  real(dp) :: R_exp(nBins)  ! Eventos medidos (Asimov = SM)
  real(dp) :: R_BSM(nBins)  ! Eventos de la nueva física
  real(dp) :: R_th(nBins)   ! Predicción teórica total (SM + BSM)
  real(dp) :: sigma2(nBins) ! Varianza estadística
  
  real(dp) :: chi2_total, alpha_min, sum_A, sum_B
  real(dp), parameter :: sigma_alpha = 0.046_dp ! 4.6% incertidumbre de flujo
  
  ! --- PARÁMETROS FIJOS ---
  real(dp), parameter :: M_S = 1.0_dp            ! Mediador fijo en 1 MeV
  real(dp), parameter :: g_S = 4.0e-6_dp         ! Acoplamiento escalar
  ! ------------------------

  ! --- PARÁMETROS DE ESCANEO PARA m_chi ---
  integer, parameter :: nSteps_mchi = 50
  real(dp) :: m_chi_val, log_mchi_min, log_mchi_max, dlog_mchi
  
  ! --- DETECTOR EFECTIVO CONUS+ ---
  real(dp), parameter :: exposure_kg_s = 1.0_dp * 119.0_dp * 86400.0_dp
  real(dp), parameter :: threshold_eV = 160.0_dp 
  real(dp), parameter :: efficiency = 1.0_dp 

  call Define_Energy_Bins(nBins, bin_edges)
  call Initialize_Flux()

  print *, "=== 1. GENERANDO DATASET DE ASIMOV (MODELO ESTANDAR) ==="
  R_exp = 0.0_dp
  
  do iBin = 1, nBins
     E_reco_min = bin_edges(iBin)
     E_reco_max = bin_edges(iBin+1)
     Eer_min_int = Eer_min_global    
     Eer_max_int = E_reco_max + 5.0_dp * sigma0
     if (Eer_min_int >= Eer_max_int) cycle

     dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
     sum_Eer = 0.0_dp

     do j = 1, nEer_points
        Eer = Eer_min_int + real(j-1, dp) * dEer
        gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
        if (gauss_int < 1.0e-10_dp) cycle

        T_N = T_from_Eer_NR(Eer)
        E_nu_min_loop = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp
        if (E_nu_min_loop > Enu_max) cycle

        dEnu = (Enu_max - E_nu_min_loop) / real(nEnu_points - 1, dp)
        sum_Enu = 0.0_dp
        
        do k = 1, nEnu_points
           E_nu = E_nu_min_loop + real(k-1, dp) * dEnu
           if (E_nu > Enu_max) E_nu = Enu_max 
           
           flux_val = Flux(E_nu)
           if (k == 1 .or. k == nEnu_points) then
              sum_Enu = sum_Enu + 0.5_dp * dSigmadEer_SM(Eer, E_nu) * flux_val
           else
              sum_Enu = sum_Enu + dSigmadEer_SM(Eer, E_nu) * flux_val
           end if
        end do
        sum_Enu = sum_Enu * dEnu * N_nuclei_per_kg * efficiency
        if (j == 1 .or. j == nEer_points) then
           sum_Eer = sum_Eer + 0.5_dp * sum_Enu * gauss_int * dEer
        else
           sum_Eer = sum_Eer + sum_Enu * gauss_int * dEer
        end if
     end do

     if ((E_reco_min * 1.0e6_dp) >= threshold_eV) then
        R_exp(iBin) = sum_Eer * exposure_kg_s
     end if
     
     sigma2(iBin) = R_exp(iBin)
     if (sigma2(iBin) == 0.0_dp) sigma2(iBin) = 1.0e-10_dp 
  end do

  print *, "=== 2. INICIANDO BARRIDO DE m_chi Y CALCULO DE CHI^2 ==="
  
  open(20, file="scan_chi2_mchi.dat", status="replace")
  write(20,'(A)') "# m_chi(MeV)       Chi^2"

  ! Escaneo logarítmico para m_chi (De 10^-3 a 10^2 MeV)
  log_mchi_min = -3.0_dp
  log_mchi_max =  2.0_dp
  dlog_mchi = (log_mchi_max - log_mchi_min) / real(nSteps_mchi - 1, dp)

  do step = 1, nSteps_mchi
     m_chi_val = 10.0_dp**(log_mchi_min + real(step-1, dp) * dlog_mchi)
     R_BSM = 0.0_dp

     ! Integración de la señal BSM para este m_chi
     do iBin = 1, nBins
        E_reco_min = bin_edges(iBin)
        E_reco_max = bin_edges(iBin+1)
        Eer_min_int = Eer_min_global    
        Eer_max_int = E_reco_max + 5.0_dp * sigma0
        if (Eer_min_int >= Eer_max_int .or. (E_reco_min * 1.0e6_dp) < threshold_eV) cycle

        dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
        sum_Eer = 0.0_dp

        do j = 1, nEer_points
           Eer = Eer_min_int + real(j-1, dp) * dEer
           gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
           if (gauss_int < 1.0e-10_dp) cycle

           ! Los umbrales cinemáticos ahora dependen dinámicamente de m_chi_val
           T_N = T_from_Eer_NR(Eer)
           Enu_min_nuc = Enu_min_Schi(T_N, m_chi_val)
           Enu_min_e = (Eer + sqrt(Eer**2 + 2.0_dp * m_e_MeV * Eer)) / 2.0_dp
           Enu_min_e = Enu_min_e * (1.0_dp + (m_chi_val**2)/(2.0_dp * m_e_MeV * Eer))

           E_nu_min_loop = min(Enu_min_nuc, Enu_min_e)
           if (E_nu_min_loop > Enu_max) cycle

           dEnu = (Enu_max - E_nu_min_loop) / real(nEnu_points - 1, dp)
           sum_Enu = 0.0_dp
           
           do k = 1, nEnu_points
              E_nu = E_nu_min_loop + real(k-1, dp) * dEnu
              if (E_nu > Enu_max) E_nu = Enu_max 
              
              flux_val = Flux(E_nu)
              if (k == 1 .or. k == nEnu_points) then
                 sum_Enu = sum_Enu + 0.5_dp * (dSigmadEer_Schi(Eer, E_nu, M_S, m_chi_val, g_S) + &
                                               dSigmadTe_Schi(Eer, E_nu, M_S, m_chi_val, g_S)) * flux_val
              else
                 sum_Enu = sum_Enu + (dSigmadEer_Schi(Eer, E_nu, M_S, m_chi_val, g_S) + &
                                      dSigmadTe_Schi(Eer, E_nu, M_S, m_chi_val, g_S)) * flux_val
              end if
           end do
           
           sum_Enu = sum_Enu * dEnu * N_nuclei_per_kg * efficiency
           if (j == 1 .or. j == nEer_points) then
              sum_Eer = sum_Eer + 0.5_dp * sum_Enu * gauss_int * dEer
           else
              sum_Eer = sum_Eer + sum_Enu * gauss_int * dEer
           end if
        end do
        R_BSM(iBin) = sum_Eer * exposure_kg_s
     end do

     ! --- CÁLCULO DEL CHI CUADRADO ---
     R_th = R_exp + R_BSM
     
     sum_A = 0.0_dp
     sum_B = 0.0_dp
     do iBin = 1, nBins
        sum_A = sum_A + ((R_th(iBin)**2) / sigma2(iBin))
        sum_B = sum_B + (R_th(iBin) * (R_exp(iBin) - R_th(iBin)) / sigma2(iBin))
     end do
     sum_A = sum_A + (1.0_dp / (sigma_alpha**2))
     alpha_min = sum_B / sum_A
     
     chi2_total = (alpha_min / sigma_alpha)**2
     do iBin = 1, nBins
        chi2_total = chi2_total + ((R_exp(iBin) - (1.0_dp + alpha_min) * R_th(iBin))**2) / sigma2(iBin)
     end do
     
     write(20,'(2ES16.5)') m_chi_val, chi2_total
     
     if (mod(step, 5) == 0 .or. step == 1) then
        print *, "Step", step, "/", nSteps_mchi, " | m_chi =", m_chi_val, "MeV | Chi^2 =", chi2_total
     end if
  end do

  close(20)
  print *, "=== ESCANEO FINALIZADO: Resultados en scan_chi2_mchi.dat ==="

end program chi2_mchi_scan
