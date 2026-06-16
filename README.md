# Metodologia_computacional_calculo_eventos_chicuadrado
# Producción de fermiones estériles en el experimento CONUS+ a partir de up-scattering

**Autor:** Esneider Isaac Guerrero Avila  
**Institución:** Universidad de Pamplona - Departamento de Física  

## Descripción del Proyecto
Este repositorio contiene las rutinas numéricas y los scripts de visualización desarrollados como parte del trabajo final de investigación para evaluar las tasas de eventos del experimento CONUS+. El marco teórico aborda la producción de fermiones estériles masivos ($\chi$) a partir de procesos inelásticos de dispersión (*up-scattering*) sobre nucleones ($\nu - N \to \chi - N$) y electrones atómicos ($\nu - e \to \chi - e$) bajo el paradigma del Portal de Neutrino Oscuro con un mediador escalar.

## Estructura del Repositorio
La arquitectura del software sigue el estándar de separación de intereses para garantizar la reproducibilidad y el orden de los cálculos:

* `src/modules/`: Contiene los módulos base en Fortran 90. Aquí se definen las constantes físicas, la cinemática, las secciones eficaces diferenciales y el procesamiento del factor de extinción de Lindhard.
* `src/main/`: Contiene los programas principales ejecutables en Fortran 90 encargados del cálculo numérico de las tasas de eventos y el análisis estadístico ($\chi^2$).
* `scripts/`: Scripts en Python orientados a la visualización de los datos resultantes, garantizando un formato de publicación profesional (con marcas internas y notación científica ajustada).

## Requisitos del Sistema
Para compilar y ejecutar los códigos de este repositorio, se requiere:
* Compilador de Fortran (ej. `gfortran`)
* Python 3.x
* Librerías de Python: `numpy`, `matplotlib`

## Instrucciones de Compilación (Fortran 90)
Dado que los programas principales dependen de los módulos de la física, la compilación debe realizarse en orden. Desde la terminal, ubicándose en la carpeta raíz del repositorio:

1. Compilar primero los módulos físicos:
   ```bash
   gfortran -c src/modules/*.f90
2. Compilar el programa principal enlazando los objetos generados
   ```bash
   gfortran src/main/nombre_del_ejecutable.f90 *.o -o calcular_tasas
3. Ejecutar el código:
    ```bash
    ./calcular_tasas
  
