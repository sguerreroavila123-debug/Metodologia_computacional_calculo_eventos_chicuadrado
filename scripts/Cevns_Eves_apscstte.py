#!/usr/bin/env python3.9
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# --- 1. Configuración de estilo "Paper" ---
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.linewidth': 2.0,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})

# --- 2. Cargar los datos ---
try:
    # 2.1 CEvNS SM
    data_cevns = np.loadtxt('/Users/hola/Documents/Fortran/Borrador/resultado.dat')
    x_bins = data_cevns[:, 0]  # Centros de bin
    y_cevns = data_cevns[:, 1]

    # 2.2 BSM Up-scattering (g_s = 2x10^-6)
    data_bsm = np.loadtxt('/Users/hola/Documents/Fortran/Borrador/espectro_CEvNS_BSM_aislado.dat')
    y_bsm = data_bsm[:, 1]

    # 2.3 BSM Up-scattering (g_s = 4x10^-6) <-- NUEVO ARCHIVO
    data_bsm_g4 = np.loadtxt('/Users/hola/Documents/Fortran/Borrador/espectro_CEvNS_BSM_aislado_gs4.dat')
    y_bsm_g4 = data_bsm_g4[:, 1]

    # 2.4 EvES SM
    data_eves = np.loadtxt('/Users/hola/Documents/Fortran/Borrador/espectro_EvES_SM.dat')
    y_eves = data_eves[:, 1]

    # 2.5 Datos experimentales CONUS+
    y_exp = np.loadtxt('data_conus_plus.dat')
    y_err_total = np.loadtxt('data_conus_plus2(1).dat')
    y_err = y_err_total / 2.0   # Dividir entre 2 para el valor ±

except FileNotFoundError as e:
    print(f"Error: No se encontró {e.filename}.")
    exit()

# --- 3. Configurar la figura ---
fig, ax = plt.subplots(figsize=(8, 6)) # Un poco más alto para hacer espacio al título

# Definir los bordes de los escalones (ancho de 10 eVee)
bin_width = 10.0
bin_edges = np.append(x_bins - bin_width/2, x_bins[-1] + bin_width/2)

# --- 4. Título de la gráfica con LaTeX ---
# Se usa r'...' para que Python interprete los comandos de LaTeX correctamente
ax.set_title(r'$M_{S} = m_{\chi} = 1$, $g_{S} = 4 \times 10^{-6}$', fontsize=16, pad=15)

# --- 5. Dibujar las contribuciones teóricas ---
# CEvNS (SM) en azul
ax.stairs(y_cevns, bin_edges, fill=False, color='royalblue', linewidth=2.0, label='CEvNS (SM)')

# BSM Up-scattering (g_s = 2x10^-6) en naranja
#ax.stairs(y_bsm, bin_edges, fill=False, color='darkorange', linewidth=2.0, label=r'BSM ($g_{S}=2\times10^{-6}$)')

# BSM Up-scattering (g_s = 4x10^-6) en morado y punteado <-- NUEVA LÍNEA
ax.stairs(y_bsm_g4, bin_edges, fill=False, color='purple', linestyle='--', linewidth=2.0, label=r'BSM ($g_{S}=4\times10^{-6}$)')

# EvES (SM) en verde
ax.stairs(y_eves, bin_edges, fill=False, color='forestgreen', linewidth=2.0, label='EvES (SM)')

# --- 6. Dibujar los datos experimentales ---
ax.errorbar(x_bins, y_exp, yerr=y_err, fmt='ko', markersize=5.5, capsize=0,
            elinewidth=1.5, label='CONUS+ data', zorder=5)

# --- 7. Añadir los umbrales de los detectores ---
thresholds = [160, 170, 180]
labels = ['C3', 'C5', 'C2']
for t, l in zip(thresholds, labels):
    ax.axvline(x=t, ymin=0, ymax=0.81, color='crimson', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.text(t, 63, l, color='crimson', ha='center', va='bottom', fontsize=11, fontfamily='serif')

# --- 8. Formateo de ejes y límites ---
ax.set_xlim(145, 350)
ax.set_ylim(-30, 80)

ax.set_xlabel(r'Energía de ionización reconstruida [eV$_{ee}$]', fontsize=16)
ax.set_ylabel(r'Excesos de eventos [kg$^{-1}$(10 eV$_{ee}$)$^{-1}$]', fontsize=16)

# Marcas mayores y menores
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1.2)

# Leyenda
ax.legend(loc='upper right', frameon=True, fontsize=11, edgecolor='black', framealpha=1)

# --- 9. Guardar y mostrar ---
plt.tight_layout()
plt.savefig('grafica_contribuciones_conus.pdf', dpi=300, bbox_inches='tight')
plt.show()
