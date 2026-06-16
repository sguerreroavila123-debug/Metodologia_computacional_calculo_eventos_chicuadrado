#!/usr/bin/env python3.9
import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar los datos del escaneo
# Asegúrate de que el archivo 'scan_chi2_MS.dat' esté en la misma carpeta
data = np.loadtxt("/Users/hola/Documents/Fortran/Borrador/scan_chi2_MS.dat")

# Extraer columnas
M_S_MeV = data[:, 0]
chi2_raw = data[:, 1]

# 2. Calcular Delta Chi^2
# Restamos el mínimo global para que la curva descanse sobre el cero
chi2_min = np.min(chi2_raw)
delta_chi2 = chi2_raw - chi2_min

# 3. Crear la gráfica
plt.figure(figsize=(8, 6))

# Graficamos la curva principal
plt.plot(M_S_MeV, delta_chi2, color='dodgerblue', linewidth=2.5, label=r'Perfil $\Delta\chi^2$')

# Añadimos la línea fundamental de exclusión al 90% CL
plt.axhline(y=2.71, color='crimson', linestyle='--', linewidth=2, label='Límite 90% CL (2.71)')

# 4. Configuración de los ejes
plt.xscale('log')
plt.ylim(0, 10)

# Etiquetas y diseño
plt.title('Análisis de Sensibilidad CONUS+\nEscaneo del Mediador Escalar ($m_\chi = 1$ MeV, $g_S = 4\\times 10^{-6}$)', fontsize=14)
plt.xlabel(r'Masa del Mediador $M_S$ [MeV]', fontsize=13)
plt.ylabel(r'$\Delta\chi^2$', fontsize=13)

# Cuadrícula fina
plt.grid(True, which="both", linestyle=':', alpha=0.7)

# --- CAMBIO REALIZADO AQUÍ ---
# Movemos la leyenda a la esquina superior izquierda
plt.legend(loc='upper left', fontsize=12)
# -----------------------------

# Ajustar bordes
plt.tight_layout()

# Mostrar la gráfica
plt.show()
