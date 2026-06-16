!================================================================================
! PROGRAMA: eventos_CEvNS_ap_aislado.f90
! Objetivo: Calcular EXCLUSIVAMENTE los eventos del nuevo fermión por up-scattering 
!           nuclear (CEvNS Escalar)
!================================================================================
program eventos_CEvNS_BSM_aislado
  use constan
  use cross_sec
  use Flujo_neu     
  use Fun_resolu
  use Integrales
  implicit none
   
  integer, parameter :: nBins = 19
  real(dp) :: bin_edges(nBins+1)
  real(dp) :: E_reco_min, E_reco_max, E_reco_center
  real(dp) :: events_BSM_bin(nBins), total_events
  integer  :: iBin, j, k
  real(dp) :: Eer, E_nu, Enu_min_req
  real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
  real(dp) :: sum_Enu, sum_Eer, gauss_int, dsigma, flux_val
  real(dp) :: exposure_kg_s, threshold_MeV, T_N
  
  ! --- PARÁMETROS BSM (Punto de prueba / Benchmark) ---
  ! Puedes modificar estos valores para ver cómo cambia el número de eventos
  real(dp), parameter :: M_S = 1.0_dp          ! Masa del mediador en MeV
  real(dp), parameter :: m_chi = 1.0_dp        ! Masa del fermión estéril en MeV
  real(dp), parameter :: g_S = 4.0e-6_dp       ! Acoplamiento escalar
  ! ----------------------------------------------------

  ! Parámetros experimentales CONUS+
  real(dp), parameter :: detector_mass_kg = 1.0_dp
  real(dp), parameter :: detector_days = 119.0_dp
  real(dp), parameter :: threshold_eV = 160.0_dp
  real(dp), parameter :: efficiency = 1.0_dp 

  call Define_Energy_Bins(nBins, bin_edges)
  threshold_MeV = threshold_eV * 1.0e-6_dp
  exposure_kg_s = detector_mass_kg * detector_days * 86400.0_dp
  call Initialize_Flux()

  print *, "=== CALCULANDO EVENTOS DE UP-SCATTERING (NUEVO FERMIÓN) ==="
  print *, "Parámetros: M_S = ", M_S, " MeV | m_chi = ", m_chi, " MeV | g_S = ", g_S
  
  open(30, file="espectro_CEvNS_BSM_aislado_gs4.dat", status="replace")
  write(30, '(A)') "# Centro_Bin(keVee)    Eventos_BSM_Puros"

  total_events = 0.0_dp

  ! BUCLE PRINCIPAL: Bines del Detector
  do iBin = 1, nBins
     E_reco_min = bin_edges(iBin)
     E_reco_max = bin_edges(iBin+1)
     E_reco_center = (E_reco_min + E_reco_max) / 2.0_dp * 1000.0_dp ! keVee

     Eer_min_int = Eer_min_global    
     Eer_max_int = E_reco_max + 5.0_dp * sigma0
     if (Eer_min_int >= Eer_max_int) cycle

     dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
     sum_Eer = 0.0_dp

     ! BUCLE DE IONIZACIÓN (Eer)
     do j = 1, nEer_points
        Eer = Eer_min_int + real(j-1, dp) * dEer
        gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
        if (gauss_int < 1.0e-10_dp) cycle

        ! Obtenemos el retroceso nuclear T_N desde la ionización medida
        T_N = T_from_Eer_NR(Eer)
        
        ! Límite cinemático estricto para crear la masa m_chi
        Enu_min_req = Enu_min_Schi(T_N, m_chi)
        if (Enu_min_req > Enu_max) cycle
        
        dEnu = (Enu_max - Enu_min_req) / real(nEnu_points - 1, dp)
        sum_Enu = 0.0_dp
        
        ! BUCLE DEL NEUTRINO (E_nu)
        do k = 1, nEnu_points
           E_nu = Enu_min_req + real(k-1, dp) * dEnu
           if (E_nu > Enu_max) E_nu = Enu_max
           flux_val = Flux(E_nu)

           ! Llamada a la sección eficaz BSM nuclear (Basada en tu ecuación de la imagen)
           ! Nota: dSigmadEer_Schi evalúa dSigmadTN_Schi y la multiplica por el Jacobiano de Quenching
           dsigma = dSigmadEer_Schi(Eer, E_nu, M_S, m_chi, g_S)

           if (k == 1 .or. k == nEnu_points) then
              sum_Enu = sum_Enu + 0.5_dp * dsigma * flux_val
           else
              sum_Enu = sum_Enu + dsigma * flux_val
           end if
        end do
        
        sum_Enu = sum_Enu * dEnu * N_nuclei_per_kg * efficiency

        if (j == 1 .or. j == nEer_points) then
           sum_Eer = sum_Eer + 0.5_dp * sum_Enu * gauss_int * dEer
        else
           sum_Eer = sum_Eer + sum_Enu * gauss_int * dEer
        end if
     end do

     events_BSM_bin(iBin) = sum_Eer * exposure_kg_s
     total_events = total_events + events_BSM_bin(iBin)
     
     write(30, '(2ES16.5)') E_reco_center, events_BSM_bin(iBin)
  end do

  print *, "--------------------------------------------------------"
  print *, "Total de Eventos del Nuevo Fermión (Up-Scattering) = ", total_events
  print *, "--------------------------------------------------------"
  
  close(30)
end program eventos_CEvNS_BSM_aislado
