#!/usr/bin/env python3.9
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar, root_scalar

# 1. Cargar y ordenar los datos 
# (CubicSpline exige que los valores de X estén en orden estrictamente creciente)
archivo = "/Users/hola/Documents/Fortran/Borrador/chi2_sin2theta.dat"
data = pd.read_csv(archivo, comment='#', sep='\s+', names=['s2w', 'chi2'])
data = data.sort_values(by='s2w')
x = data['s2w'].values
y = data['chi2'].values

# 2. Interpolación Cúbica de los datos brutos
# Esto crea una función continua que pasa EXACTAMENTE por tus puntos de Fortran
spline = CubicSpline(x, y)

# 3. Encontrar el mínimo EXACTO numéricamente
# Buscamos el mínimo de la función interpolada dentro del rango de tus datos
res = minimize_scalar(spline, bounds=(x.min(), x.max()), method='bounded')
min_s2w = res.x
y_min = res.fun

print(f"--- Resultados del Análisis con Interpolación ---")
print(f"Mejor valor sin²θ_W (mínimo real): {min_s2w:.6f}")
print(f"Chi² mínimo original: {y_min:.6f}")

# 4. Definir función Δχ² restando el mínimo
def delta_chi2(val):
    return spline(val) - y_min

# 5. Calcular los límites de 1-sigma (Δχ² = 1)
# Creamos una función auxiliar que sea 0 cuando Δχ² = 1
def func_root(val):
    return delta_chi2(val) - 1.0

# Buscamos la raíz izquierda (x < min_s2w)
root_left = root_scalar(func_root, bracket=[x.min(), min_s2w]).root
# Buscamos la raíz derecha (x > min_s2w)
root_right = root_scalar(func_root, bracket=[min_s2w, x.max()]).root

print(f"Intervalo 1-sigma: [{root_left:.6f}, {root_right:.6f}]")
error_izq = min_s2w - root_left
error_der = root_right - min_s2w
print(f"Incertidumbre: -{error_izq:.6f} / +{error_der:.6f}")

# 6. Visualización
plt.figure(figsize=(8, 6))

# Crear una malla fina para que la línea dibujada sea perfectamente suave
x_smooth = np.linspace(x.min(), x.max(), 1000)
y_shifted_smooth = delta_chi2(x_smooth)

# Graficamos la curva interpolada ya restada
plt.plot(x_smooth, y_shifted_smooth, '-', label='Interpolación Δχ²', color='blue', linewidth=2)

# Opcional: Superponemos los puntos reales de tu archivo para que tu tutora vea
# que la curva pasa exactamente por los puntos calculados en tu simulación
plt.plot(x, y - y_min, 'ko', markersize=3, label='Datos Fortran')

plt.axhline(1.0, color='red', linestyle='--', label=r'Δ$\chi^2 = 1$ (1$\sigma$)')
plt.axvline(min_s2w, color='green', linestyle=':', label=f'Mínimo: {min_s2w:.3f}')

# Restricciones de la gráfica (Petición de la tutora)
plt.ylim(0, 10)  
plt.xlim(min_s2w - 0.15, min_s2w + 0.15) 

plt.xlabel(r'$\sin^2\theta_W$')
plt.ylabel(r'$\Delta\chi^2$') 
plt.title(r'Análisis de sensibilidad de $\sin^2\theta_W$ (CONUS+)')
plt.legend()
plt.grid(True)
plt.show()
