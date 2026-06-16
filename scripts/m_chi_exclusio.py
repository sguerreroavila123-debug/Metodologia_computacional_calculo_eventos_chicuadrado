#!/usr/bin/env python3.9
import numpy as np
import matplotlib.pyplot as plt

# 1. Configuración de estilo "Paper"
plt.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.linewidth': 1.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})

# 2. Función para cargar y estructurar la matriz 2D
def procesar_matriz(filename, n_pts=100):
    data = np.loadtxt(filename)
    mchi_flat = data[:, 0]
    gS_flat = data[:, 1]
    chi2_flat = data[:, 2]

    # Convertir a grilla 2D
    mchi_grid = mchi_flat.reshape(n_pts, n_pts)
    gS_grid = gS_flat.reshape(n_pts, n_pts)
    chi2_grid = chi2_flat.reshape(n_pts, n_pts)

    # Calcular Delta Chi^2 restando el fondo (mínimo absoluto)
    chi2_min = np.min(chi2_grid)
    delta_chi2_grid = chi2_grid - chi2_min
    
    return mchi_grid, gS_grid, delta_chi2_grid

# 3. Cargar los dos escenarios
# Asegúrate de que los archivos .dat estén en el mismo directorio
mchi_nuc, gS_nuc, chi2_nuc = procesar_matriz("/Users/hola/Documents/Fortran/Borrador/matriz_chi2_2D_mchi.dat")
mchi_comp, gS_comp, chi2_comp = procesar_matriz("/Users/hola/Documents/Fortran/Borrador/matriz_chi2_2D_mchi_Completa.dat")

# 4. Crear la gráfica
fig, ax = plt.subplots(figsize=(9, 6))

# Trazar el escenario Completo (Línea Roja y Sombreado)
ax.contour(mchi_comp, gS_comp, chi2_comp, levels=[2.71], 
           colors='red', linewidths=2.5)
ax.contourf(mchi_comp, gS_comp, chi2_comp, levels=[2.71, 1e10], 
            colors=['red'], alpha=0.1)

# Trazar el escenario Solo Nuclear (Línea Azul Punteada)
ax.contour(mchi_nuc, gS_nuc, chi2_nuc, levels=[2.71], 
           colors='blue', linestyles='--', linewidths=2.5)

# Configuración de los ejes (Escala Logarítmica)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 1e2)
ax.set_ylim(1e-7, 1e-2)

# Etiquetas
ax.set_xlabel(r'Masa del fermión estéril $m_\chi$ [MeV]', fontsize=16)
ax.set_ylabel(r'Acoplamiento escalar $g_S$', fontsize=16)
#ax.set_title(r'Comparativa de Sensibilidad ($M_S = 1.0$ MeV)', fontsize=16)

# Leyenda manual
ax.plot([], [], color='red', linewidth=2.5, 
        label=r'CE$\nu$NS (SM) + CE$\nu$NS($\chi$) + EvES ($\chi$)')
ax.plot([], [], color='blue', linestyle='--', linewidth=2.5, 
        label=r'CE$\nu$NS (SM) + CE$\nu$NS ($\chi$)')
ax.legend(loc='upper left', fontsize=12, framealpha=1.0, edgecolor='black')

# Estética
ax.grid(True, which="both", linestyle=':', alpha=0.6)

plt.tight_layout()
plt.show()
