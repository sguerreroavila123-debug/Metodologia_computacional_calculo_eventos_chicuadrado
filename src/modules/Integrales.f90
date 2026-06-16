!================================================================================
! MÓDULO: Integrales_M
! Define los bines de energía reconstruida (E_reco) para comparar con los
! excesos de conteos de CONUS+ (19 bines, 160 - 350 eVee).
!================================================================================
module Integrales
  use constan
  implicit none
  
contains

  !----------------------------------------------------------------------------
  ! Genera los límites de los bines de energía en MeV
  ! n = 19 para la replicación exacta de CONUS+
  !----------------------------------------------------------------------------
  subroutine Define_Energy_Bins(n, bin_edges)
    integer, intent(in) :: n
    real(dp), intent(out) :: bin_edges(n+1)
    integer :: i
    real(dp) :: E_start_eV, bin_width_eV
    
    E_start_eV = 160.0_dp
    bin_width_eV = 10.0_dp
    
    do i = 1, n+1
       ! Casteo explícito de 'i' para evitar pérdida de precisión
       bin_edges(i) = (E_start_eV + real(i-1, dp) * bin_width_eV) * 1.0e-6_dp   ! eVee -> MeV
    end do
    
    print *, "========================================="
    print *, "BINS DE ENERGÍA RECONSTRUIDA (E_reco)"
    print *, "========================================="
    do i = 1, n
       print '(A,I2,A,F6.1,A,F6.1,A)', " Bin ", i, ": ", bin_edges(i)*1.0e6_dp, &
                                       " - ", bin_edges(i+1)*1.0e6_dp, " eVee"
    end do
    print *, "========================================="
  end subroutine Define_Energy_Bins

  !----------------------------------------------------------------------------
  ! Calcula el centro exacto de un bin (útil para gráficas)
  !----------------------------------------------------------------------------
  pure real(dp) function Bin_Center(bin_edges, i)
    real(dp), intent(in) :: bin_edges(:)
    integer, intent(in) :: i
    
    Bin_Center = (bin_edges(i) + bin_edges(i+1)) / 2.0_dp
  end function Bin_Center

end module Integrales
