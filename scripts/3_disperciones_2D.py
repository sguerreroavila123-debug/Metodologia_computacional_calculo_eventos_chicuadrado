OA#!/usr/bin/env python3.9
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# --- 1. Configuración de estilo "Paper" ---
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.linewidth': 1.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})

# --- 2. Cargar los eventos del Modelo Estándar ---
sm_events = np.array([33.04, 25.67, 19.97, 15.55, 12.10, 9.40, 7.30, 5.65, 4.36, 
    3.35, 2.57, 1.96, 1.48, 1.12, 0.84, 0.63, 0.47, 0.35, 0.26])

# Asumimos que la varianza es igual al número de eventos (estadística de Poisson)
sigma2_sm = sm_events.copy() 

# --- 3. FUNCIÓN PARA CALCULAR LÍMITES DE EXCLUSIÓN ---
# Esta función hace el trabajo duro: lee el archivo y calcula la curva
def calcular_limites_exclusion(archivo_matriz, g_ref=1.0e-4):
    try:
        bsm_data = np.loadtxt(archivo_matriz)
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo {archivo_matriz}")
        return None, None
        
    M_S_vals = bsm_data[:, 0]
    bsm_spectra = bsm_data[:, 2:]
    
    limites_gS = [] 
    
    for i in range(len(M_S_vals)):
        N_BSM_ref = bsm_spectra[i]
        
        # Si no hay eventos BSM, no hay sensibilidad
        if np.sum(N_BSM_ref) <= 0:
            limites_gS.append(np.nan)
            continue

        # Definir la función Delta Chi^2 para buscar su raíz (Delta Chi^2 = 2.71 para 90% CL)
        def delta_chi2(g_S):
            factor_escala = (g_S / g_ref)**4
            N_BSM_test = N_BSM_ref * factor_escala
            chi2 = np.sum((N_BSM_test**2) / sigma2_sm)
            return chi2 - 2.71

        try:
            sol = root_scalar(delta_chi2, bracket=[1e-8, 15.0], method='brentq')
            limites_gS.append(sol.root)
        except ValueError:
            limites_gS.append(np.nan)
            
    return M_S_vals, np.array(limites_gS)

# --- 4. PROCESAR LOS DOS ARCHIVOS ---
# Archivo 1
archivo_1 = "/Users/hola/Documents/Fortran/Borrador/matriz_BSM_espectro33.dat"
MS_1, limites_1 = calcular_limites_exclusion(archivo_1)

# Archivo 2 (Reemplaza con tu ruta si es diferente)
archivo_2 = "/Users/hola/Documents/Fortran/Borrador/matriz_BSM_Completa3escenarios.dat"
MS_2, limites_2 = calcular_limites_exclusion(archivo_2)

# --- 5. GRAFICAR ---
fig, ax = plt.subplots(figsize=(8, 6))

# Dibujar el Escenario 1 (Azul)
if MS_1 is not None:
    ax.plot(MS_1, limites_1, color='blue', linewidth=2.5, label=r'CE$\nu$NS (SM) + CE$\nu$NS ($\chi$)')
    ax.fill_between(MS_1, limites_1, 1e-2, color='royalblue', alpha=0.15)

# Dibujar el Escenario 2 (Rojo/Carmesí)
if MS_2 is not None:
    # Usamos linestyle='--' para diferenciar mejor si las líneas se cruzan
    ax.plot(MS_2, limites_2, color='crimson', linewidth=2.5, linestyle='--', label=r'CE$\nu$NS (SM) + CE$\nu$NS ($\chi$) + EvES ($\chi$)')
    # Usamos un alpha distinto para que al superponerse se vean ambos colores
    ax.fill_between(MS_2, limites_2, 1e-2, color='lightcoral', alpha=0.15)

# Configuración de ejes logarítmicos y límites
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 1e3)
ax.set_ylim(1e-7, 1e-2)

# Textos en español
ax.set_xlabel(r'Masa del mediador $M_S$ [MeV]', fontsize=16)
ax.set_ylabel(r'Acoplamiento escalar $g_S$', fontsize=16)
#ax.set_title(r'Proyección CONUS+: $\it{Up-scattering}$ de fermión estéril', fontsize=14, pad=15)

# Cuadrícula y leyenda
ax.grid(True, which="major", linestyle="-", alpha=0.4)
ax.grid(True, which="minor", linestyle=":", alpha=0.4)
ax.legend(loc='upper left', frameon=True, fontsize=11, edgecolor='black')

# Guardar y mostrar
plt.tight_layout()
plt.savefig('grafica_exclusion_comparativa.pdf', dpi=300, bbox_inches='tight')
plt.show()
