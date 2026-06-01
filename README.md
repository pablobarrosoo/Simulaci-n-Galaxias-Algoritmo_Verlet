# Simulacion-Galaxias-Algoritmo_Verlet
Proyecto de Física Computacional para explicar la forma espiral de las galaxias

## Programas en C++:
Los programas del método 1 están listos para ejecutar y observar la simulación. Si el programa de cálculos tarda demasiado recomiendo bajar el número de cuerpos N. Esto se hace en la cabecera de todos los programas, en la parte de Datos.

Los programas del método 2 tiene dos configuraciones. Actualmente se encuentra configurado para observar la evolución de los brazos espirales. En el caso de que se quiera usar el programa para el cálculo del periodo del Sol deben hacerse algunos cambios en el apartado de Datos al inicio del código:
1) Disminuir considerablemente el número de cuerpos. Recomiendo 500 o 100.
2) Aumentar las iteraciones a 4000.
3) Elegir el semi eje mayor del Sol deseado. Cota superior o inferior.
4) Seleccionar la excentricidad del Sol igual a 0,06, la real.


## Programas de Python:
Se encargan de representar las simulaciones. Está todo listo para ejecutar. Solo hay que asegurarse de que el número de cuerpos (N) sea el mismo que en C++. Para acelerar o ralentizar la simulación recomiendo modificar el parametro de 'SALTO_FRAMES', 1 para que dibuje todos los frames y X para que dibuje 1 de cada x frames.
