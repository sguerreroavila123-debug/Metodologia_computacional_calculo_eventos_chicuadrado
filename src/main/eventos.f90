!================================================================================
! PROGRAMA PRINCIPAL:  (con flujo combinado Kopeikin+Mueller)
! Calcula la predicción SM de eventos CEvNS para el detector efectivo de CONUS+
! (1 kg, 119 días, umbral 160 eVee)
!================================================================================
program eventos
  use constan
  use cross_sec
  use Flujo_neu     
  use Fun_resolu
  use Integrales
  implicit none
   
  integer, parameter :: nBins = 19
  real(dp) :: bin_edges(nBins+1)
  real(dp) :: E_reco_min, E_reco_max
  real(dp) :: events_bin(nBins), total_events
  integer  :: iBin, j, k
  real(dp) :: Eer, E_nu, E_nu_min, T_N
  real(dp) :: Eer_min_int, Eer_max_int, dEer, dEnu
  real(dp) :: sum_Enu, sum_Eer, gauss_int, dsigma, flux_val
  real(dp) :: exposure_kg_s, threshold_MeV
  
  ! Parámetros experimentales específicos para el detector
  real(dp), parameter :: detector_mass_kg = 1.0_dp
  real(dp), parameter :: detector_days = 119.0_dp
  real(dp), parameter :: threshold_eV = 160.0_dp
  real(dp), parameter :: efficiency = 1.0_dp 

  call Define_Energy_Bins(nBins, bin_edges)

  threshold_MeV = threshold_eV * 1.0e-6_dp
  exposure_kg_s = detector_mass_kg * detector_days * 86400.0_dp
  
  ! Inicializar flujo combinado (Asumiendo que Flujo_neu expone esta subrutina)
  call Initialize_Flux()

  print *, ""
  print *, "=== INICIANDO INTEGRACIÓN NUMÉRICA CEνNS (SM)  ==="
  
  total_events = 0.0_dp
  events_bin = 0.0_dp

  ! Bucle sobre los bines de energía reconstruida (E_reco)
  do iBin = 1, nBins
     E_reco_min = bin_edges(iBin)
     E_reco_max = bin_edges(iBin+1)

     Eer_min_int = Eer_min_global    
     Eer_max_int = E_reco_max + 5.0_dp * sigma0
     if (Eer_min_int >= Eer_max_int) cycle

     dEer = (Eer_max_int - Eer_min_int) / real(nEer_points - 1, dp)
     sum_Eer = 0.0_dp

     ! Bucle sobre la energía verdadera de ionización (Eer)
     do j = 1, nEer_points
        Eer = Eer_min_int + real(j-1, dp) * dEer
        
        gauss_int = Integral_Gaussiana(E_reco_min, E_reco_max, Eer)
        if (gauss_int < 1.0e-10_dp) cycle

        ! Límite cinemático EXACTO para no evaluar zonas de sección eficaz nula
        T_N = T_from_Eer_NR(Eer)
        E_nu_min = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp
        if (E_nu_min > Enu_max) cycle

        dEnu = (Enu_max - E_nu_min) / real(nEnu_points - 1, dp)
        sum_Enu = 0.0_dp
        
        ! Bucle sobre la energía del neutrino (E_nu)
        do k = 1, nEnu_points
           E_nu = E_nu_min + real(k-1, dp) * dEnu
           
           ! Eliminamos el 'exit' de punto flotante para no truncar la regla del trapecio
           if (E_nu > Enu_max) E_nu = Enu_max 

           dsigma = dSigmadEer_SM(Eer, E_nu)
           flux_val = Flux(E_nu)

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

     events_bin(iBin) = sum_Eer * exposure_kg_s
     total_events = total_events + events_bin(iBin)

     write(*,'(A,I2,A,F8.2,A,ES10.3)') " Bin ", iBin, " (centro ", &
                Bin_Center(bin_edges, iBin)*1.0e6_dp, " eV): ", events_bin(iBin)
  end do

  print *, "====================================================="
  print *, "Total eventos SM :", total_events
  print *, "====================================================="

  open(10, file="resultado.dat")
  do iBin = 1, nBins
     write(10,'(F8.2,2X,ES15.5)') Bin_Center(bin_edges, iBin)*1.0e6_dp, events_bin(iBin)
  end do
  close(10)
  print *, "Resultados guardados en resultado.dat"
end program eventos
