#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <chrono>
#include <random>
#include <ctime>

using namespace std;




///////////////////////////////////////////////////////////////////////////////////////////////////////

//Datos
#define M_sagitarioA 8.546733e36            //Masa de Sagitario A* en kg.

//Definimos constantes de reescalado
#define Ms 1.9891e30                        //Masa del Sol. Será la unidad de masa en el sistema reescalado
#define f_x (1.0/3.0857e19)                 //Factor de escala de la distancia. Unidades de distancia: Kiloparsec (kpc).
#define f_t (1.0/3.1536e13)                 //Factor de escala del tiempo. Unidades de tiempo: Millones de años (Myr).
#define f_v (1.022004732e-6)                //Factor de escala de la velocidad. Unidades de velocidad: kpc/Myr.

//Constantes para la simulación
#define N 2000                              //Número de cuerpos a simular.
#define iteraciones_max 1000                //Número máximo de iteraciones.
#define salto_t 1e1                         //Salto de tiempo. 100000años. 1e-1Myr.

//Constantes reescaladas
#define G_reescalada 4.493379481e-12        //Constante de gravitación reescalada, calculada a partir de las constantes de reescalado y la constante de gravitación original
#define M_sagitarioA_reescalada 4.297e6     //Masa de Sagitario A* reescalada, calculada a partir de la masa original y el factor de escala de la masa.

//Constantes para la espiral
#define a_espiral 2                         //Constante a de la espiral logarítmica
#define b_espiral 0.3                       //Constante b de la espiral logarítmica
#define Amp_ruido 0.5                       //Amplitud del ruido aleatorio que se le añade a la posición de cada cuerpo para simular la distribución de los cuerpos en una galaxia real
#define pi 3.14159265358979323846           //Valor de pi




///////////////////////////////////////////////////////////////////////////////////////////////////////

//DECLARACION DE FUNCIONES

///////////////////////////////////////////////////////////////////////////////////////////////////////

//Funciones para rellenar matrices y vectores
void rellena_m (double m[]);                                                //Función para rellenar el vector m. Todas tendra 1 masa solar
void rellena_inicial_r (double r[][2], mt19937_64& generator);              //Función para rellenar la matriz r con las posiciones iniciales de los cuerpos, distribuidos en una espiral logaritmica con ruido aleatorio
void rellena_v(double r[][2], double v[][2]);                               //Función para rellenar la matriz v con las velocidades iniciales, que se calculan a partir de las posiciones iniciales y la masa de Sagitario A*


//Funciones para trabajar con matrices y vectores
void rellena_cero_matriz (double a[][2]);                                   //Inicializa a cero una matriz


//Funciones para realizar operaciones matemáticas
void resta_vec (double r[][2], double r_resta[], int v1, int v2);           //Función para restar dos vectores
double mod_vec (double v[]);                                                //Función para calcular modulo de un vector


//Funciones para realizar el algoritmo
void calculo_a (double r[][2], double a[][2], double m[]);                  //Función para calcular a
void calculo_r (double r[][2], double v[][2], double a[][2], double h);     //Funcion para calcular r
void calculo_w (double w[][2], double v[][2], double a[][2], double h);     //Funcion para calcular w
void calculo_v (double w[][2], double v[][2], double a[][2], double h);     //Funcion para calcular v


//Funciones para números aleatorios
double r_azar (mt19937_64& generator);                                      //Función para generar número aleatorio real entre 0 y 1




///////////////////////////////////////////////////////////////////////////////////////////////////////

//MAIN

///////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){

    /////////////////////////////////////////////////////////////////////////
    //Definimos variables

    double m[N], r[N][2], v[N][2], a[N][2], w[N][2];

    double t = 0;    //Tiempo inicial




    /////////////////////////////////////////////////////////////////////////
    //Generamos la semilla para los numeros aleatorios

    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);




    /////////////////////////////////////////////////////////////////////////
    //Rellenamos las matrices y vectores con los datos iniciales

    //Vector de masas:
    m[0] = M_sagitarioA_reescalada;
    rellena_m(m);


    //Vector de posiciones iniciales:
    r[0][0] = 0;   //Posición inicial de Sagitario A* (cuerpo 0)
    r[0][1] = 0;
    rellena_inicial_r(r, generator);


    //Vector de velocidades iniciales:
    v[0][0] = 0;   //Velocidad inicial de Sagitario A* (cuerpo 0)
    v[0][1] = 0;
    rellena_v(r, v);




    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //Abirmos fichero para guardar resultados
    ofstream arch_resultados("resultados.txt"); //Abrimos archivo de resultados
    ofstream arch_resultados_simp("resultados_simp.txt");   //Abrimos archivo de resultados simplificado

    if(!arch_resultados){   //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de resultados" << endl;
    }

    if(!arch_resultados_simp){  //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de resultados" << endl;
    }

    //Encabezado de los archivos
    arch_resultados << "Tiempo" << "\t" << "Cuerpo" << "\t" << "rx" << "\t" << "ry" << "\t" << "vx" << "\t" << "vy" << "\t" << "ax" << "\t" << "ay" << endl;
    arch_resultados_simp << "Tiempo" << "\t" << "Cuerpo" << "\t" << "rx" << "\t" << "ry" << endl;



    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //Algoritmo

    //Inicializamos la aceleración a cero y caluclamos la aceleración por primera vez, el resto de veces será en el bucle
    rellena_cero_matriz(a);
    calculo_a(r,a,m);

    int n=0;
    while (n<iteraciones_max){

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //Imprimimos resultados antes de sobrescribir

        for (int i=0; i<N; i++){
            arch_resultados << t << "\t" << i << "\t";
            arch_resultados_simp << t << "\t"<< i << "\t";
            
            for (int j=0; j<2; j++){
                arch_resultados << r[i][j] << "\t";
                arch_resultados_simp << r[i][j] << "\t";
            }

            for (int j=0; j<2; j++){
                arch_resultados << v[i][j] << "\t";
            }

            for (int j=0; j<2; j++){
                arch_resultados << a[i][j] << "\t";
            }

            arch_resultados << endl;
            arch_resultados_simp << endl;
        }




        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //Procedimiento del algoritmo
        calculo_r(r,v,a,salto_t);       //Calculo de r(t+h)
        calculo_w(w,v,a,salto_t);       //Calculo de w, donde se guarda valores de v(t) y a(t)

        rellena_cero_matriz(a);         //Inicializamos a
        calculo_a (r,a,m);              //Calculamos a(t+h)
        calculo_v(w,v,a,salto_t);       //Calculamos v(t+h)




        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //Para conteo, guardamos el tiempo y contamos iteraciones
        cout << "Iteracion: " << n << "\r";
        t += salto_t;
        n += 1;
    }
    
    return 0;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////

//DEFINIMOS FUNCIONES

///////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para rellenar matrices y vectores

//Función para rellenar el vector masa
void rellena_m (double m[]){
    for (int i=1; i<N; i++){
        m[i] = 1;
    }
}


//Función para rellenar la matriz r con las posiciones iniciales de los cuerpos, distribuidos en una espiral logaritmica con ruido aleatorio
void rellena_inicial_r(double r[][2], mt19937_64& generator){

    double theta, r_espiral, ruido, r_final;

    //Brazo 1
    for (int i=1; i<N/2; i++){

        theta = 2*pi*r_azar(generator);                     //Ángulo aleatorio entre 0 y 2pi

        r_espiral = a_espiral*exp(b_espiral*theta);         //Radio de la espiral logarítmica para el ángulo theta

        ruido = 2*Amp_ruido*r_azar(generator) - Amp_ruido;  //Ruido aleatorio entre -Amp_ruido y Amp_ruido

        r_final = r_espiral + ruido;                        //Radio final, que es la suma del radio de la espiral y el ruido aleatorio

        r[i][0] = r_final*cos(theta);       //Coordenada x
        r[i][1] = -r_final*sin(theta);      //Coordenada y
    }

    //Brazo 2
    for (int i=N/2; i<N; i++){

        theta = 2*pi*r_azar(generator);                     //Ángulo aleatorio entre 0 y 2pi

        r_espiral = a_espiral*exp(b_espiral*theta);         //Radio de la espiral logarítmica para el ángulo theta

        ruido = 2*Amp_ruido*r_azar(generator) - Amp_ruido;  //Ruido aleatorio entre -Amp_ruido y Amp_ruido

        r_final = r_espiral + ruido;                        //Radio final, que es la suma del radio de la espiral y el ruido aleatorio

        r[i][0] = r_final*cos(theta + pi);      //Coordenada x
        r[i][1] = -r_final*sin(theta + pi);     //Coordenada y
    }
}


//Función para rellenar la matriz v con las velocidades iniciales, que se calculan a partir de las posiciones iniciales y la masa de Sagitario A*
void rellena_v(double r[][2], double v[][2]){

    double modulo_r, v_mag;

    for (int i=1; i<N; i++){

        modulo_r = mod_vec(r[i]);   //Calculamos el módulo de la posición del cuerpo i

        v_mag = sqrt(G_reescalada*M_sagitarioA_reescalada/modulo_r);   //Calculamos la velocidad tangencial necesaria para que el cuerpo i esté en órbita circular alrededor de Sagitario A* a la distancia dada por su posición inicial.

        v[i][0] = -v_mag*r[i][1]/modulo_r;   //Calculamos la componente x de la velocidad, que es perpendicular a la posición y tiene magnitud v_mag
        v[i][1] = v_mag*r[i][0]/modulo_r;    //Calculamos la componente y de la velocidad, que es perpendicular a la posición y tiene magnitud v_mag
    }
}


//Funcion para inicializar a cero
void rellena_cero_matriz (double a[][2]){
    for (int i=0; i<N; i++){
        for (int j=0; j<2; j++){
            a[i][j] = 0;
        }
    }
}




///////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para operaciones matemáticas

//Restamos v1 menos v2
void resta_vec (double r[][2], double r_resta[], int v1, int v2){ 
    for (int j=0; j<2; j++){
        r_resta[j] = r[v1][j] - r[v2][j];
    }
}

//Módulo de un vector
double mod_vec (double v[]){
    double mod;
    
    mod = sqrt(v[0]*v[0] + v[1]*v[1]);

    return mod;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para el algoritmo

//Todas las operaciones empiezan a partir del cuerpo 1, que es el que orbita alrededor de Sagitario A* (cuerpo 0). El cuerpo 0 se mantiene fijo en el origen, por lo que no se calcula su movimiento.

//Calculo de la aceleración
void calculo_a (double r[][2], double a[][2], double m[]){

    double r_resta[2];
    double mod; 

    for (int i=1; i<N; i++){

        //Algunas operaciones auxiliares
        resta_vec(r, r_resta, i,0);
        mod = mod_vec(r_resta);

        for (int k=0; k<2; k++) {   //Cada componente de la aceleración
            a[i][k] += -G_reescalada*m[0]*r_resta[k]/pow(mod,3);
        }
    }   
}



//Calculo de la posición
void calculo_r (double r[][2], double v[][2], double a[][2], double h){
    for (int i=1; i<N; i++){
            for (int j=0; j<2; j++){
                r[i][j] = r[i][j] + h*v[i][j] + h*h*0.5*a[i][j];
            }
        }
}



//Calculo del vector auxiliar w
void calculo_w (double w[][2], double v[][2], double a[][2], double h){
    for (int i=1; i<N; i++){
        for (int j=0; j<2; j++){
            w[i][j] = v[i][j] + h*0.5*a[i][j]; //Vector auxiliar que guarda a(t) y que luego será usado en el calculo de v(t+h)
        }
    }
}



//Calculo de la velocidad
void calculo_v (double w[][2], double v[][2], double a[][2], double h){
    for (int i=1; i<N; i++){
            for (int j=0; j<2; j++){
                v[i][j] = w[i][j] + h*0.5*a[i][j];
            }
        }
}




///////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para numero aleatorios

//Función para generar número aleatorio real entre 0 y 1
double r_azar (mt19937_64& generator) {

    uniform_real_distribution<double> r_distribution(0., 1.); //initialize the distribution r_distribution
    
    return r_distribution(generator); //a random double extracted from the distribution r_distribution
}

