!================================================================================
! PROGRAMA: chi2_2D_mchi_gs_completo.f90
! Objetivo: Escaneo 2D de Chi^2 en el espacio de parámetros (m_chi, g_S)
! Física: CEνNS (SM) + CEνNS (BSM) + EvES (BSM). CANAL ELECTRÓNICO ACTIVADO.
! Condición: M_S fijo (ej. 1.0 MeV), variando m_chi y g_S.
!================================================================================
program chi2_2D_mchi_gs_completo
  use constan
  use cross_sec
  use Flujo_neu     
  use Fun_resolu
  use Integrales
  implicit none
   
  integer, parameter :: nBins = 19
  real(dp) :: bin_edges(nBins+1)
  real(dp) :: E_reco_min, E_reco_max
  
  ! Arreglos de eventos
  real(dp) :: R_exp(nBins), R_BSM_ref(nBins), R_th(nBins), sigma2(nBins)
  real(dp) :: R_BSM_scaled(nBins)
  
  integer  :: iBin, j, k, step_m, step_g
  real(dp) :: Eer, E_nu, T_N
  real(dp) :: Enu_min_SM, Enu_min_BSM_nuc, Enu_min_e, Enu_min_BSM_e, E_nu_min_loop
  real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
  real(dp) :: sum_Enu_SM, sum_Enu_BSM, sum_Eer_SM, sum_Eer_BSM
  real(dp) :: gauss_int, dsigma_SM, dsigma_BSM_nuc, dsigma_BSM_e, dsigma_BSM, flux_val
  
  ! Variables Chi Cuadrado
  real(dp) :: chi2_total, alpha_min, sum_A, sum_B
  real(dp), parameter :: sigma_alpha = 0.046_dp ! 4.6% de incertidumbre
  
  ! --- PARÁMETROS DEL ESCANEO 2D ---
  integer, parameter :: nSteps_mchi = 100
  integer, parameter :: nSteps_gS = 100
  
  real(dp) :: m_chi_val, log_mchi_min, log_mchi_max, dlog_mchi
  real(dp) :: gS_val, log_gS_min, log_gS_max, dlog_gS
  
  ! --- PARÁMETRO FIJO ---
  real(dp), parameter :: MS_val = 1.0_dp        ! M_S fijo en 1.0 MeV
  real(dp), parameter :: g_ref = 1.0e-4_dp      ! Acoplamiento base para integrar rápido
  real(dp) :: factor_escala
  
  ! Detector efectivo CONUS+
  real(dp), parameter :: exposure_kg_s = 1.0_dp * 119.0_dp * 86400.0_dp
  real(dp), parameter :: threshold_eV = 160.0_dp 
  real(dp), parameter :: efficiency = 1.0_dp 

  call Define_Energy_Bins(nBins, bin_edges)
  call Initialize_Flux()

  ! ==================================================================
  ! PASO 1: Calcular Fondo Modelo Estándar (R_exp)
  ! ==================================================================
  print *, "=== 1. CALCULANDO FONDO DEL MODELO ESTANDAR (CEvNS) ==="
  R_exp = 0.0_dp
  
  do iBin = 1, nBins
     E_reco_min = bin_edges(iBin)
     E_reco_max = bin_edges(iBin+1)
     Eer_min_int = Eer_min_global    
     Eer_max_int = E_reco_max + 5.0_dp * sigma0
     if (Eer_min_int >= Eer_max_int) cycle

     dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
     sum_Eer_SM = 0.0_dp

     do j = 1, nEer_points
        Eer = Eer_min_int + real(j-1, dp) * dEer
        gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
        if (gauss_int < 1.0e-10_dp) cycle

        T_N = T_from_Eer_NR(Eer)
        Enu_min_SM = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp
        if (Enu_min_SM > Enu_max) cycle

        dEnu = (Enu_max - Enu_min_SM) / real(nEnu_points - 1, dp)
        sum_Enu_SM = 0.0_dp
        
        do k = 1, nEnu_points
           E_nu = Enu_min_SM + real(k-1, dp) * dEnu
           if (E_nu > Enu_max) E_nu = Enu_max 
           flux_val = Flux(E_nu)
           
           if (k == 1 .or. k == nEnu_points) then
              sum_Enu_SM = sum_Enu_SM + 0.5_dp * dSigmadEer_SM(Eer, E_nu) * flux_val
           else
              sum_Enu_SM = sum_Enu_SM + dSigmadEer_SM(Eer, E_nu) * flux_val
           end if
        end do
        sum_Enu_SM = sum_Enu_SM * dEnu * N_nuclei_per_kg * efficiency
        if (j == 1 .or. j == nEer_points) then
           sum_Eer_SM = sum_Eer_SM + 0.5_dp * sum_Enu_SM * gauss_int * dEer
        else
           sum_Eer_SM = sum_Eer_SM + sum_Enu_SM * gauss_int * dEer
        end if
     end do
     
     if ((E_reco_min * 1.0e6_dp) >= threshold_eV) then
        R_exp(iBin) = sum_Eer_SM * exposure_kg_s
     end if
     sigma2(iBin) = R_exp(iBin)
     if (sigma2(iBin) == 0.0_dp) sigma2(iBin) = 1.0e-10_dp 
  end do

  ! ==================================================================
  ! PASO 2: BUCLE 2D SOBRE MASA DEL FERMIÓN (m_chi) Y ACOPLAMIENTOS (g_S)
  ! ==================================================================
  print *, "=== 2. INICIANDO ESCANEO 2D: (m_chi, g_S) (Nuclear + Electrónico) ==="
  
  ! Cambiamos el nombre del archivo para no sobreescribir el anterior
  open(20, file="matriz_chi2_2D_mchi_Completa.dat", status="replace")
  write(20,'(A)') "# m_chi(MeV)       g_S             Chi2"

  log_mchi_min = -3.0_dp
  log_mchi_max =  2.0_dp
  dlog_mchi = (log_mchi_max - log_mchi_min) / real(nSteps_mchi - 1, dp)

  log_gS_min = -7.0_dp
  log_gS_max = -2.0_dp
  dlog_gS = (log_gS_max - log_gS_min) / real(nSteps_gS - 1, dp)

  do step_m = 1, nSteps_mchi
     m_chi_val = 10.0_dp**(log_mchi_min + real(step_m-1, dp) * dlog_mchi)
     R_BSM_ref = 0.0_dp

     ! Filtro cinemático conservador basado en la energía del reactor
     if (m_chi_val < Enu_max) then
         
         do iBin = 1, nBins
            E_reco_min = bin_edges(iBin)
            E_reco_max = bin_edges(iBin+1)
            Eer_min_int = Eer_min_global    
            Eer_max_int = E_reco_max + 5.0_dp * sigma0
            if (Eer_min_int >= Eer_max_int .or. (E_reco_min * 1.0e6_dp) < threshold_eV) cycle

            dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
            sum_Eer_BSM = 0.0_dp

            do j = 1, nEer_points
               Eer = Eer_min_int + real(j-1, dp) * dEer
               gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
               if (gauss_int < 1.0e-10_dp) cycle

               ! 1. Umbral Estéril NUCLEAR
               T_N = T_from_Eer_NR(Eer)
               Enu_min_SM = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp
               Enu_min_BSM_nuc = Enu_min_SM * (1.0_dp + (m_chi_val**2)/(2.0_dp * M_nucleus_Ge * T_N))
               
               ! 2. Umbral Estéril ELECTRÓNICO
               Enu_min_e = (Eer + sqrt(Eer**2 + 2.0_dp * m_e_MeV * Eer)) / 2.0_dp
               Enu_min_BSM_e = Enu_min_e * (1.0_dp + (m_chi_val**2)/(2.0_dp * m_e_MeV * Eer))
               
               ! El flujo debe empezar desde el canal que se active primero
               E_nu_min_loop = min(Enu_min_BSM_nuc, Enu_min_BSM_e)
               if (E_nu_min_loop > Enu_max) cycle

               dEnu = (Enu_max - E_nu_min_loop) / real(nEnu_points - 1, dp)
               sum_Enu_BSM = 0.0_dp
               
               do k = 1, nEnu_points
                  E_nu = E_nu_min_loop + real(k-1, dp) * dEnu
                  if (E_nu > Enu_max) E_nu = Enu_max 
                  flux_val = Flux(E_nu)
                  
                  dsigma_BSM_nuc = 0.0_dp
                  dsigma_BSM_e = 0.0_dp
                  
                  ! Evaluar canal nuclear si E_nu es suficiente
                  if (E_nu >= Enu_min_BSM_nuc) then
                     dsigma_BSM_nuc = dSigmadEer_Schi(Eer, E_nu, MS_val, m_chi_val, g_ref)
                  end if
                  
                  ! Evaluar canal electrónico si E_nu es suficiente
                  if (E_nu >= Enu_min_BSM_e) then
                     dsigma_BSM_e = dSigmadTe_Schi(Eer, E_nu, MS_val, m_chi_val, g_ref)
                  end if
                  
                  ! Suma total de secciones eficaces BSM
                  dsigma_BSM = dsigma_BSM_nuc + dsigma_BSM_e

                  if (k == 1 .or. k == nEnu_points) then
                     sum_Enu_BSM = sum_Enu_BSM + 0.5_dp * dsigma_BSM * flux_val
                  else
                     sum_Enu_BSM = sum_Enu_BSM + dsigma_BSM * flux_val
                  end if
               end do
               sum_Enu_BSM = sum_Enu_BSM * dEnu * N_nuclei_per_kg * efficiency
               if (j == 1 .or. j == nEer_points) then
                  sum_Eer_BSM = sum_Eer_BSM + 0.5_dp * sum_Enu_BSM * gauss_int * dEer
               else
                  sum_Eer_BSM = sum_Eer_BSM + sum_Enu_BSM * gauss_int * dEer
               end if
            end do
            R_BSM_ref(iBin) = sum_Eer_BSM * exposure_kg_s
         end do
     end if

     ! Bucle interno: escalar g_S y calcular Chi^2
     do step_g = 1, nSteps_gS
        gS_val = 10.0_dp**(log_gS_min + real(step_g-1, dp) * dlog_gS)
        
        ! Ambos canales escalan con g_S^4
        factor_escala = (gS_val / g_ref)**4
        R_BSM_scaled = R_BSM_ref * factor_escala
        R_th = R_exp + R_BSM_scaled
        
        ! Minimización
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
        
        ! Escribir (m_chi, g_S, chi2)
        write(20,'(3ES16.5)') m_chi_val, gS_val, chi2_total
     end do
     
     if (mod(step_m, 10) == 0) print *, "Progreso m_chi: ", step_m, "/", nSteps_mchi
  end do

  close(20)
  print *, "=== FINALIZADO. Datos exportados a matriz_chi2_2D_mchi_Completa.dat ==="
end program chi2_2D_mchi_gs_completo
