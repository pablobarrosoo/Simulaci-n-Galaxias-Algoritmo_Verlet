import numpy as np
import matplotlib.pyplot as plt

############################################
# PARÁMETROS DEL MODELO
############################################
N = 500             # Número de elipses (órbitas de estrellas)
ecc = 0.5           # Excentricidad fija (según tu enunciado)
a_min = 1.0         # Semieje mayor mínimo (kpc)
a_max = 12.0        # Semieje mayor máximo (kpc)
limite = 15         # Límite de los ejes para la gráfica (kpc)

# Factor de torsión: cuánto se rota el eje de la elipse por cada unidad de semieje mayor (en radianes)
k_torsion = -0.6        

# Configurar el estilo visual
plt.figure(figsize=(9, 9))
plt.style.use('dark_background')

############################################
# GENERACIÓN DE DATOS
############################################
# Generamos semiejes mayores 'a' aleatorios entre el rango mínimo y máximo
a_valores = np.random.uniform(a_min, a_max, N)


for i, a in enumerate(a_valores):

    # Elegimos el ángulo de rotación en función del semieje mayor de las elipses
    if (i<N/2):
        phi = k_torsion * a  - k_torsion * a_min
    else: 
        phi = k_torsion * a  - k_torsion * a_min + np.pi
    
    # Generamos el ángulo pola
    theta = np.linspace(0, 2 * np.pi, 300)
    
    # Ecuación de la elipse con foco en el origen
    r = (a * (1 - ecc**2)) / (1 - ecc * np.cos(theta - phi))
    
    # Transformación de coordenadas polares a cartesianas (X, Y)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    # Dibujamos la línea de la órbita
    plt.plot(x, y, color='yellow', alpha=0.15, linewidth=1.2)
    
    

############################################
# CONFIGURACIÓN DEL GRÁFICO
############################################
# Dibujamos el agujero negro
plt.scatter(0, 0, color='white', s=40, zorder=5, label='Sagitario A*')

# Ajustes de los ejes
plt.xlim(-limite, limite)
plt.ylim(-limite, limite)
plt.gca().set_aspect('equal', adjustable='box')

plt.title("Representación de elipses rotadas", fontsize=25, pad=15, color='black')
plt.xlabel("Eje X (kpc)", color='black', fontsize=20)
plt.ylabel("Eje Y (kpc)", color='black', fontsize=20)
plt.tick_params(axis='both', colors='black', labelsize=15, bottom=True, left=True, labelbottom=True, labelleft=True)
plt.grid(color='gray', linestyle='--', alpha=0.1)
plt.legend(loc='upper right')

plt.show()