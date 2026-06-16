!================================================================================
! PROGRAMA PRINCIPAL
! FÍSICA:analisi \chi^{2} 2d (Ms,gs) CEνNS (SM) + Up-Scattering Nuclear (BSM) + EvES (BSM)
!================================================================================
program eventos_BSM_Completo
  use constan
  use cross_sec
  use Flujo_neu     
  use Fun_resolu
  use Integrales
  implicit none
   
  integer, parameter :: nBins = 19
  real(dp) :: bin_edges(nBins+1)
  real(dp) :: E_reco_min, E_reco_max
  real(dp) :: events_SM_bin(nBins), events_BSM_bin(nBins), events_Total_bin(nBins)
  real(dp) :: total_events_BSM
  integer  :: iBin, j, k, iMass
real(dp) :: Eer, E_nu, E_nu_min_SM, E_nu_min_BSM_nuc, E_nu_min_BSM_e, E_nu_min_loop
  real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
  real(dp) :: sum_Enu_SM, sum_Enu_BSM, sum_Eer_SM, sum_Eer_BSM
  real(dp) :: gauss_int, dsigma_SM, dsigma_BSM_nuc, dsigma_BSM_e, flux_val
  real(dp) :: exposure_kg_s, threshold_MeV, T_N
  
  ! --- PARÁMETROS BSM ---
  integer, parameter :: nMassPoints = 200      
  real(dp) :: M_S, m_chi, m_chi_max
  real(dp) :: log_MS_min, log_MS_max, dlog_MS
  real(dp), parameter :: g_ref = 1.0d-4        
  real(dp), parameter :: ratio_mass = 1.0_dp   
  ! ----------------------------------------------

  ! Parámetros experimentales C3
  real(dp), parameter :: detector_mass_kg = 1.0_dp
  real(dp), parameter :: detector_days = 119.0_dp
  real(dp), parameter :: threshold_eV = 160.0_dp
  real(dp), parameter :: efficiency = 1.0_dp 

  call Define_Energy_Bins(nBins, bin_edges)
  threshold_MeV = threshold_eV * 1.0e-6_dp
  exposure_kg_s = detector_mass_kg * detector_days * 86400.0_dp
  call Initialize_Flux()

  print *, "=== BARRIDO: CEνNS (SM) + UP-SCATTERING (NUCLEAR + ELECTRÓNICO) ==="
  
  open(20, file="matriz_BSM_Completa3escenarios.dat", status="replace")
  
  log_MS_min = -3.0_dp
  log_MS_max =  3.0_dp
  dlog_MS = (log_MS_max - log_MS_min) / real(nMassPoints - 1, dp)

  do iMass = 1, nMassPoints
     M_S = 10.0_dp**(log_MS_min + real(iMass - 1, dp) * dlog_MS)
     m_chi = ratio_mass * M_S
     
     ! Límite cinemático simple
     m_chi_max = sqrt(M_nucleus_Ge * (M_nucleus_Ge + 2.0_dp * Enu_max)) - M_nucleus_Ge
     if (m_chi >= m_chi_max) cycle
     
     events_SM_bin = 0.0_dp
     events_BSM_bin = 0.0_dp
     total_events_BSM = 0.0_dp

     do iBin = 1, nBins
        E_reco_min = bin_edges(iBin)
        E_reco_max = bin_edges(iBin+1)
        Eer_min_int = Eer_min_global    
        Eer_max_int = E_reco_max + 5.0_dp * sigma0
        if (Eer_min_int >= Eer_max_int) cycle

        dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
        sum_Eer_SM = 0.0_dp
        sum_Eer_BSM = 0.0_dp

        do j = 1, nEer_points
           Eer = Eer_min_int + real(j-1, dp) * dEer
           gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
           if (gauss_int < 1.0e-10_dp) cycle

           T_N = T_from_Eer_NR(Eer)
           
           ! Cinemática de los 3 canales
           E_nu_min_SM = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp
           E_nu_min_BSM_nuc = E_nu_min_SM * (1.0_dp + (m_chi**2)/(2.0_dp * M_nucleus_Ge * T_N))
           
           ! Umbral EvES para la energía de retroceso actual
           E_nu_min_BSM_e = (Eer + sqrt(Eer**2 + 2.0_dp * m_e_MeV * Eer)) / 2.0_dp
           E_nu_min_BSM_e = E_nu_min_BSM_e * (1.0_dp + (m_chi**2)/(2.0_dp * m_e_MeV * Eer))

           ! Definimos el inicio de la integración como el mínimo de los canales activos
           E_nu_min_loop = min(E_nu_min_SM, E_nu_min_BSM_nuc, E_nu_min_BSM_e)
           if (E_nu_min_loop > Enu_max) cycle
           
           dEnu = (Enu_max - E_nu_min_loop) / real(nEnu_points - 1, dp)
           sum_Enu_SM = 0.0_dp
           sum_Enu_BSM = 0.0_dp
           
           do k = 1, nEnu_points
              E_nu = E_nu_min_loop + real(k-1, dp) * dEnu
              if (E_nu > Enu_max) E_nu = Enu_max
     
              flux_val = Flux(E_nu)

              ! 1. Sección Eficaz SM (Solo contribuye si E_nu > SM threshold)
              dsigma_SM = 0.0_dp
              if (E_nu >= E_nu_min_SM) dsigma_SM = dSigmadEer_SM(Eer, E_nu)

              ! 2. Sección Eficaz BSM Nuclear (CEvNS)
              dsigma_BSM_nuc = 0.0_dp
              if (E_nu >= E_nu_min_BSM_nuc) then
                 dsigma_BSM_nuc = dSigmadEer_Schi(Eer, E_nu, M_S, m_chi, g_ref)
              end if

              ! 3. Sección Eficaz BSM Electrónico (EvES) - AHORA ACTIVADO
              dsigma_BSM_e = 0.0_dp
              if (E_nu >= E_nu_min_BSM_e) then
                 dsigma_BSM_e = dSigmadTe_Schi(Eer, E_nu, M_S, m_chi, g_ref)
              end if

              ! Integración sumando todos los términos BSM
              if (k == 1 .or. k == nEnu_points) then
                 sum_Enu_SM  = sum_Enu_SM  + 0.5_dp * dsigma_SM * flux_val
                 sum_Enu_BSM = sum_Enu_BSM + 0.5_dp * (dsigma_BSM_nuc + dsigma_BSM_e) * flux_val
              else
                 sum_Enu_SM  = sum_Enu_SM  + dsigma_SM * flux_val
                 sum_Enu_BSM = sum_Enu_BSM + (dsigma_BSM_nuc + dsigma_BSM_e) * flux_val
              end if
           end do
           
           sum_Enu_SM  = sum_Enu_SM * dEnu * N_nuclei_per_kg * efficiency
           sum_Enu_BSM = sum_Enu_BSM * dEnu * N_nuclei_per_kg * efficiency

           if (j == 1 .or. j == nEer_points) then
              sum_Eer_SM  = sum_Eer_SM  + 0.5_dp * sum_Enu_SM * gauss_int * dEer
              sum_Eer_BSM = sum_Eer_BSM + 0.5_dp * sum_Enu_BSM * gauss_int * dEer
           else
              sum_Eer_SM  = sum_Eer_SM  + sum_Enu_SM * gauss_int * dEer
              sum_Eer_BSM = sum_Eer_BSM + sum_Enu_BSM * gauss_int * dEer
           end if
        end do

        events_SM_bin(iBin) = sum_Eer_SM * exposure_kg_s
        events_BSM_bin(iBin) = sum_Eer_BSM * exposure_kg_s
        total_events_BSM = total_events_BSM + events_BSM_bin(iBin)
     end do

     print *, "M_S = ", M_S, " | Eventos BSM Totales = ", total_events_BSM
     write(20, '(2ES14.5, 19ES14.5)') M_S, m_chi, events_BSM_bin
  end do

  close(20)
  print *, "=== BARRIDO COMPLETADO (CANALES NUCLEAR + ELECTRÓNICO) ==="
end program eventos_BSM_Completo
