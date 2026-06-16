!================================================================================
! MÓDULO: cross_sec
!   Contiene: Lindhard QF, derivada analítica, resolución T_N desde Eer (NR),
!   factor de forma de Klein-Nystrand, sección eficaz SM (CEνNS) y EνES.
!   Totalmente migrado al estándar paramétrico real(dp).
!   CORREGIDO: Consistencia dimensional en la derivada de Lindhard (Jacobiano).
!================================================================================
module cross_sec
  use constan
  implicit none

contains

  !----------------------------------------------------------------------------
  ! Factor de Forma Nuclear de Klein-Nystrand (Paper Ec. 3)
  ! Recibe T_N en MeV
  !----------------------------------------------------------------------------
  pure real(dp) function FormFactor_KN2(T_N)
    real(dp), intent(in) :: T_N
    real(dp) :: q_MeV, q_fm, R_A, a_k, x, j1
    
    ! Si la energía es ~ 0, el factor es 1
    if (T_N <= 1.0e-10_dp) then
       FormFactor_KN2 = 1.0_dp
       return
    end if

    ! Transferencia de momento q = sqrt(2 * M * T_N) (en MeV)
    q_MeV = sqrt(2.0_dp * M_nucleus_Ge * T_N)
    ! Convertir momento a fm^-1 dividiendo por hbarc_fm (~197.32 MeV fm)
    q_fm = q_MeV / hbarc_fm
    
    R_A = 1.23_dp * (A_Ge**(1.0_dp/3.0_dp)) ! fm
    a_k = 0.7_dp ! fm
    x = q_fm * R_A
    
    ! Función esférica de Bessel de orden 1
    j1 = (sin(x) / (x**2)) - (cos(x) / x)
    
    ! Factor de forma al cuadrado
    FormFactor_KN2 = ( 3.0_dp * j1 / x * (1.0_dp / (1.0_dp + (q_fm**2 * a_k**2))) )**2
  end function FormFactor_KN2

  !----------------------------------------------------------------------------
  ! Lindhard Quenching Factor (recibe T_N en MeV)
  !----------------------------------------------------------------------------
  pure real(dp) function Lindhard_QF(T_N_MeV)
    real(dp), intent(in) :: T_N_MeV
    real(dp) :: T_N_keV, eps, g_eps

    if (T_N_MeV <= 0.0_dp) then
       Lindhard_QF = 0.2_dp
       return
    end if

    T_N_keV = T_N_MeV * 1000.0_dp
    eps = 11.5_dp * T_N_keV * (Z_Lindhard)**(-7.0_dp/3.0_dp)
    g_eps = 3.0_dp * (eps**0.15_dp) + 0.7_dp * (eps**0.6_dp) + eps
    Lindhard_QF = (k_Lindhard * g_eps) / (1.0_dp + k_Lindhard * g_eps)
  end function Lindhard_QF

  !----------------------------------------------------------------------------
  ! Derivada analítica de QF respecto a T_N (en MeV⁻¹) - ¡CORREGIDO!
  !----------------------------------------------------------------------------
  pure real(dp) function dQF_dT_analytic(T_N_MeV)
    real(dp), intent(in) :: T_N_MeV
    real(dp) :: T_keV, eps, g, dg_deps, alpha_keV
    
    T_keV = T_N_MeV * 1000.0_dp
    
    ! La constante empírica alfa evalúa T en keV
    alpha_keV = 11.5_dp * (Z_Lindhard)**(-7.0_dp/3.0_dp)
    eps = alpha_keV * T_keV
    
    g = 3.0_dp * (eps**0.15_dp) + 0.7_dp * (eps**0.6_dp) + eps
    dg_deps = 0.45_dp * (eps**(-0.85_dp)) + 0.42_dp * (eps**(-0.4_dp)) + 1.0_dp
    
    ! Calculamos la derivada respecto a T_keV y multiplicamos por 1000 
    ! al final para obtener la dimensión correcta en MeV⁻¹
    dQF_dT_analytic = (k_Lindhard * dg_deps * alpha_keV / (1.0_dp + k_Lindhard * g)**2) * 1000.0_dp
  end function dQF_dT_analytic

  !----------------------------------------------------------------------------
  ! Resuelve T_N a partir de Eer usando Newton-Raphson (más rápido y robusto)
  !----------------------------------------------------------------------------
  pure real(dp) function T_from_Eer_NR(Eer)
    real(dp), intent(in) :: Eer
    real(dp) :: T, f, df, delta
    integer :: iter
    integer, parameter :: max_iter = 50
    real(dp), parameter :: tol = 1.0e-12_dp

    if (Eer <= 0.0_dp) then
       T_from_Eer_NR = 0.0_dp
       return
    end if

    ! Estimación inicial inteligente basada en el rango de CONUS+
    if (Eer < 1.0e-4_dp) then
       T = Eer * 4.5_dp       ! QF bajo a muy bajas energías
    else
       T = Eer * 5.0_dp       ! QF medio a altas energías
    end if

    do iter = 1, max_iter
       f  = T * Lindhard_QF(T) - Eer
       df = Lindhard_QF(T) + T * dQF_dT_analytic(T)

       if (abs(df) < 1.0e-14_dp) exit   ! Evita división por cero severa

       delta = f / df
       T = T - delta
       
       ! Restricción física (evitar retrocesos negativos)
       if (T < 0.0_dp) T = Eer * 5.0_dp ! Reset en caso de oscilación salvaje

       if (abs(delta) < tol) exit
    end do

    T_from_Eer_NR = T
  end function T_from_Eer_NR

  !----------------------------------------------------------------------------
  ! dσ/dT_N del SM para CEνNS (Paper Ec. 1)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadTN_SM(T_N, E_nu)
    real(dp), intent(in) :: T_N, E_nu
    real(dp) :: Q_V, gVp, gVn, F_W_sq

    ! Restricción Cinemática
    if (E_nu < sqrt(M_nucleus_Ge * T_N) / 2.0_dp) then
       dSigmadTN_SM = 0.0_dp
       return
    end if

    gVp = 0.5_dp - 2.0_dp * sin2W
    gVn = -0.5_dp
    Q_V = gVp * Z_Ge + gVn * N_Ge
    F_W_sq = FormFactor_KN2(T_N)

    dSigmadTN_SM = (GF_MeV**2 * M_nucleus_Ge / pi) * (Q_V**2) * 1.0_dp * &
                   (1.0_dp - ((M_nucleus_Ge * T_N) / (2.0_dp * E_nu**2)) - (T_N / E_nu))
                   
    if (dSigmadTN_SM < 0.0_dp) dSigmadTN_SM = 0.0_dp
  end function dSigmadTN_SM

  !----------------------------------------------------------------------------
  ! dσ/dEer CEνNS usando la ecuación 28 del paper (con derivada analítica de QF)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadEer_SM(Eer, E_nu)
    real(dp), intent(in) :: Eer, E_nu
    real(dp) :: T_N, QF, dQF_dT_val, dT_dEer, dsigma_dTN

    T_N = T_from_Eer_NR(Eer)
    if (T_N <= 0.0_dp) then
       dSigmadEer_SM = 0.0_dp
       return
    end if
    
    QF = Lindhard_QF(T_N)
    dQF_dT_val = dQF_dT_analytic(T_N)
    
    ! Jacobiano dT_N/dEer = 1 / (dEer/dT_N) donde Eer = T_N * QF(T_N)
    dT_dEer = 1.0_dp / (QF + T_N * dQF_dT_val)
    dsigma_dTN = dSigmadTN_SM(T_N, E_nu)
    
    ! Se aplica conv_MeV2_to_cm2 para que el resultado esté en cm² / MeV
    dSigmadEer_SM = dsigma_dTN * dT_dEer * conv_MeV2_to_cm2
  end function dSigmadEer_SM

  !----------------------------------------------------------------------------
  ! Sección Eficaz dσ/dEer del SM para Dispersión de Electrones (EνES) (Paper Ec. 4)
  !----------------------------------------------------------------------------

  !

  
  

pure real(dp) function Enu_min_Schi(T_N, m_chi)
    real(dp), intent(in) :: T_N       ! Retroceso nuclear en MeV
    real(dp), intent(in) :: m_chi     ! Masa del fermión estéril en MeV
    real(dp) :: Enu_min_SM, factor_masa

    if (T_N <= 0.0_dp) then
       Enu_min_Schi = Enu_max ! Si no hay retroceso, requiere energía infinita
       return
    end if

    ! 1. Calculamos el límite estándar del SM (m_chi = 0) que ya conoces
    Enu_min_SM = (T_N + sqrt(T_N**2 + 2.0_dp * M_nucleus_Ge * T_N)) / 2.0_dp

    ! 2. Aplicamos la penalización por la creación de la masa estéril
    factor_masa = 1.0_dp + (m_chi**2) / (2.0_dp * M_nucleus_Ge * T_N)

    Enu_min_Schi = factor_masa * Enu_min_SM
  end function Enu_min_Schi
  
  !----------------------------------------------------------------------------
  ! BSM: Masa máxima permitida del fermión estéril por límite del reactor
  ! (Basado en la sección H del artículo de CONUS+)
  !----------------------------------------------------------------------------
  pure real(dp) function max_m_chi_allowed(E_nu_max_local)
    real(dp), intent(in) :: E_nu_max_local
    max_m_chi_allowed = sqrt(M_nucleus_Ge * (M_nucleus_Ge + 2.0_dp * E_nu_max_local)) - M_nucleus_Ge
  end function max_m_chi_allowed

  !----------------------------------------------------------------------------
  ! BSM: dσ/dT_N para CEνNS con Mediador Escalar y Fermión Estéril (Ec. 22)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadTN_Schi(T_N, E_nu, m_S, m_chi, g_S)
    real(dp), intent(in) :: T_N, E_nu, m_S, m_chi, g_S
    real(dp) :: factor_p, factor_n, Q_S_nuc, F_W_sq, Enu_min_req, propagador

    ! 1. Verificación cinemática exacta (Fase 1)
    Enu_min_req = Enu_min_Schi(T_N, m_chi)
    if (E_nu < Enu_min_req) then
       dSigmadTN_Schi = 0.0_dp
       return
    end if

    ! 2. Carga Escalar Nuclear Q_S (Basado en Ec. 11 con g_S universal)
    factor_p = (m_p_MeV / m_u_MeV)*f_Tu_p + (m_p_MeV / m_d_MeV)*f_Td_p
    factor_n = (m_n_MeV / m_u_MeV)*f_Tu_n + (m_n_MeV / m_d_MeV)*f_Td_n
    Q_S_nuc = (g_S**2) * (Z_Ge * factor_p + N_Ge * factor_n)

    ! 3. Factor de forma de Klein-Nystrand
    F_W_sq = FormFactor_KN2(T_N)

    ! 4. Propagador masivo del mediador escalar
    propagador = 1.0_dp / (m_S**2 + 2.0_dp * M_nucleus_Ge * T_N)**2

    ! 5. Ensamblaje de la Ecuación 22
    dSigmadTN_Schi = (M_nucleus_Ge * Q_S_nuc**2) / (4.0_dp * pi) * propagador * F_W_sq &
                   * (1.0_dp + T_N / (2.0_dp * M_nucleus_Ge)) &
                   * ( (M_nucleus_Ge * T_N)/(E_nu**2) + (m_chi**2)/(2.0_dp * E_nu**2) )
                   
    if (dSigmadTN_Schi < 0.0_dp) dSigmadTN_Schi = 0.0_dp
  end function dSigmadTN_Schi

  !----------------------------------------------------------------------------
  ! BSM: Transformación a dσ/dEer para CEνNS Escalar (Incluye Quenching)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadEer_Schi(Eer, E_nu, m_S, m_chi, g_S)
    real(dp), intent(in) :: Eer, E_nu, m_S, m_chi, g_S
    real(dp) :: T_N, QF_val, dQF_dT_val, dT_dEer, dsigma_dTN

    T_N = T_from_Eer_NR(Eer)
    if (T_N <= 0.0_dp) then
       dSigmadEer_Schi = 0.0_dp
       return
    end if
    
    QF_val = Lindhard_QF(T_N)
    dQF_dT_val = dQF_dT_analytic(T_N)
    
    dT_dEer = 1.0_dp / (QF_val + T_N * dQF_dT_val)
    dsigma_dTN = dSigmadTN_Schi(T_N, E_nu, m_S, m_chi, g_S)
    
    ! Conversión final a cm²/MeV
    dSigmadEer_Schi = dsigma_dTN * dT_dEer * conv_MeV2_to_cm2
  end function dSigmadEer_Schi

  !----------------------------------------------------------------------------
  ! BSM: dσ/dEer para Dispersión de Electrones (EνES) con Mediador Escalar (Ec. 24)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadTe_Schi(Eer, E_nu, m_S, m_chi, g_S)
    real(dp), intent(in) :: Eer, E_nu, m_S, m_chi, g_S
    real(dp) :: Enu_min_req, Z_eff, propagador

    ! 1. Cinemática para el electrón (debe tener energía para crear m_chi)
    Enu_min_req = (Eer + sqrt(Eer**2 + 2.0_dp * m_e_MeV * Eer)) / 2.0_dp
    Enu_min_req = Enu_min_req * (1.0_dp + (m_chi**2)/(2.0_dp * m_e_MeV * Eer))
    
    if (E_nu < Enu_min_req) then
       dSigmadTe_Schi = 0.0_dp
       return
    end if

    ! Aproximación de Z_eff a baja energía (citado en el paper)
    Z_eff = Z_Ge

    ! 2. Propagador masivo con la masa del electrón
    propagador = 1.0_dp / (m_S**2 + 2.0_dp * m_e_MeV * Eer)**2

    ! 3. Ensamblaje de la Ecuación 24
    dSigmadTe_Schi = (m_e_MeV * g_S**4) / (4.0_dp * pi) * propagador * Z_eff &
                   * (1.0_dp + Eer / (2.0_dp * m_e_MeV)) &
                   * ( (m_e_MeV * Eer)/(E_nu**2) + (m_chi**2)/(2.0_dp * E_nu**2) )
                   
    ! Conversión final a cm²/MeV
    dSigmadTe_Schi = dSigmadTe_Schi * conv_MeV2_to_cm2
  end function dSigmadTe_Schi

 !----------------------------------------------------------------------------
  ! Sección Eficaz dσ/dEer del SM para Dispersión de Electrones (EνES)
  !----------------------------------------------------------------------------
  pure real(dp) function dSigmadTe_SM(T, E_nu)
    real(dp), intent(in) :: T, E_nu
    real(dp) :: gV, gA, dsig

    ! Acoplamientos efectivos para antineutrinos electrónicos (Interferencia CC + NC)
    gV = 2.0_dp * 0.23867_dp + 0.5_dp
    gA = -0.5_dp

    ! Usando GF_MeV y conv_MeV2_to_cm2 de tu módulo constan.f90
    dsig = 32.0_dp * (GF_MeV**2) * m_e_MeV / (2.0_dp * pi) * &
           ( (gV+gA)**2 + (gV-gA)**2 * (1.0_dp - T/E_nu)**2 - (gV**2 - gA**2) * (m_e_MeV*T)/(E_nu**2) )

    ! Convertir unidades
    dSigmadTe_SM = dsig * conv_MeV2_to_cm2
  end function dSigmadTe_SM


  
end module cross_sec
