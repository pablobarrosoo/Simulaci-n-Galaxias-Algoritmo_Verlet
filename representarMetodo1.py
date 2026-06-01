import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#######################################################################################################

#Datos de entrada y constantes
ARCHIVO_DATOS = "resultados_simp.txt"       # Nombre del archivo exportado por tu programa en C++
N = 2000                                    # Numero de cuerpos (IMPORTANTE: Debe ser el mismo N que en C++)

#Constantes para la grafica y animacion
LIMITE_EJES = 17.5                          # Limite de distancia (kpc)
INTERVALO_FRAME = 20                        # Milisegundos de pausa entre cada iteracion
SALTO_FRAMES =  10                          # Dibuja solo 1 de cada X iteraciones




#######################################################################################################
print("Cargando...")

#Leemos los datos del archivo generado por C++
   
# Saltamos la primera linea (encabezado) y leemos las 4 columnas
datos = np.loadtxt(ARCHIVO_DATOS, skiprows=1)   # Estructura: [Tiempo, Cuerpo, rx, ry]

# Extraemos los tiempos unicos para saber cuantas iteraciones totales tenemos
tiempos = np.unique(datos[:, 0])
num_frames_totales = len(tiempos)
    
# Aplicamos un salto para no procesar el 100% de los frames y que la animacion sea fluida
indices_frames = np.arange(0, num_frames_totales, SALTO_FRAMES)
tiempos_anim = tiempos[indices_frames]
num_frames = len(indices_frames)

# Reestructuramos los datos en matrices para acceder facilmente por indice [frame, cuerpo]
x_matrix = datos[:, 2].reshape(num_frames_totales, N)
y_matrix = datos[:, 3].reshape(num_frames_totales, N)




#########################################################################
#Configuramos la figura y los ejes de la ventana

fig, ax = plt.subplots(figsize=(8, 8))        
    
# Formateamos el aspecto estético del grafico
ax.set_facecolor('black')
ax.set_xlim(-LIMITE_EJES, LIMITE_EJES)
ax.set_ylim(-LIMITE_EJES, LIMITE_EJES)
ax.set_xlabel("Eje X (kpc)", fontsize='20')
ax.set_ylabel("Eje Y (kpc)", fontsize='20')
ax.set_title("Problema del 'Winding'", fontsize='25')
ax.tick_params(axis='both', labelsize=12)



#########################################################################
#Inicializamos los elementos visuales

# Sagitario A* (cuerpo 0)
agujero_negro, = ax.plot([], [], 'o', color='white', markersize=5, label="Sagitario A*")
    
# Estrellas (cuerpos del 1 al N-1)
estrellas, = ax.plot([], [], '.', color='yellow', markersize=4, alpha=0.7, label="Estrellas")

ax.legend(loc="upper right")

# Cuadro de texto para mostrar el tiempo en pantalla
texto_tiempo = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)




#########################################################################
#Definimos las funciones del bucle de la animacion

# Funcion para limpiar el grafico al inicio
def inicializa_grafico():
    agujero_negro.set_data([], [])
    estrellas.set_data([], [])
    texto_tiempo.set_text('')
    return agujero_negro, estrellas, texto_tiempo

# Funcion que actualiza las posiciones en cada frame de la animacion
def actualiza(frame_idx):
    # Obtenemos el indice real del tiempo considerando los saltos de frames
    n = indices_frames[frame_idx]

    # Posicion de Sagitario A* (cuerpo indice 0)
    x_sagA = [x_matrix[n, 0]]
    y_sagA = [y_matrix[n, 0]]
    agujero_negro.set_data(x_sagA, y_sagA)

    # Posiciones del resto de estrellas (cuerpos indices del 1 al N-1)
    x_estr = x_matrix[n, 1:]
    y_estr = y_matrix[n, 1:]
    estrellas.set_data(x_estr, y_estr)

    # Actualizamos el rotulo del tiempo actual
    t_actual = tiempos_anim[frame_idx]
    texto_tiempo.set_text(f"Tiempo: {t_actual:.2f} Myr")

    return agujero_negro, estrellas, texto_tiempo




#########################################################################
#Ejecutamos la animacion
    
anim = animation.FuncAnimation(fig, actualiza, frames=num_frames, init_func=inicializa_grafico, blit=True, interval=INTERVALO_FRAME)

# Mostramos la animacion por pantalla
plt.show()