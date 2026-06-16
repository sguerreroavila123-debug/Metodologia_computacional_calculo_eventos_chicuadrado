#!/usr/bin/env python3.9

import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar los datos del nuevo escaneo
# Asegúrate de que el archivo 'scan_chi2_mchi.dat' esté en la misma carpeta
data = np.loadtxt("/Users/hola/Documents/Fortran/Borrador/scan_chi2_mchi.dat")

# Extraer columnas
m_chi_MeV = data[:, 0]
chi2_raw = data[:, 1]

# 2. Calcular Delta Chi^2
# Restamos el mínimo global para anclar la curva al cero
chi2_min = np.min(chi2_raw)
delta_chi2 = chi2_raw - chi2_min

# 3. Crear la gráfica
plt.figure(figsize=(8, 6))

# Graficamos la curva principal
plt.plot(m_chi_MeV, delta_chi2, color='darkorange', linewidth=2.5, label=r'Perfil $\Delta\chi^2$')

# Añadimos la línea fundamental de exclusión al 90% CL
plt.axhline(y=2.71, color='crimson', linestyle='--', linewidth=2, label='Límite 90% CL (2.71)')

# 4. Configuración de los ejes
plt.xscale('log')
plt.ylim(0, 10)

# Etiquetas y diseño actualizados para m_chi
plt.title('Análisis de Sensibilidad CONUS+\nEscaneo de la Masa del Fermión ($M_S = 1$ MeV, $g_S = 4\\times 10^{-6}$)', fontsize=14)
plt.xlabel(r'Masa del Fermión Estéril $m_\chi$ [MeV]', fontsize=13)
plt.ylabel(r'$\Delta\chi^2$', fontsize=13)

# Cuadrícula fina
plt.grid(True, which="both", linestyle=':', alpha=0.7)

# Leyenda en la esquina superior izquierda
plt.legend(loc='upper left', fontsize=12)

# Ajustar bordes
plt.tight_layout()

# Mostrar la gráfica
plt.show()
