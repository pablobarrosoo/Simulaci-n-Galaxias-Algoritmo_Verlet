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
#define M_sagitarioA 8.546733e36   //Masa de Sagitario A* en kg.

//Definimos constantes de reescalado
#define Ms 1.9891e30                //Masa del Sol. Será la unidad de masa en el sistema reescalado
#define f_x (1.0/3.0857e19)         //Factor de escala de la distancia. Unidades de distancia: Kiloparsec (kpc).
#define f_t (1.0/3.1536e13)         //Factor de escala del tiempo. Unidades de tiempo: Millones de años (Myr).
#define f_v (1.022004732e-6)        //Factor de escala de la velocidad. Unidades de velocidad: kpc/Myr.

//Constantes para la simulación
#define N 5000                      //Número de cuerpos a simular.
#define iteraciones_max 1000        //Número máximo de iteraciones.
#define salto_t 1e1                 //Salto de tiempo. 100000años.

//Constantes reescaladas
#define G_reescalada 4.493379481e-12        //Constante de gravitación reescalada, calculada a partir de las constantes de reescalado y la constante de gravitación original
#define M_sagitarioA_reescalada 4.297e6     //Masa de Sagitario A* reescalada, calculada a partir de la masa original y el factor de escala de la masa.

//Constantes para las elipses
#define a_min 1                     //Semieje mayor mínimo de las elipses en kpc
#define a_max 12                    //Semieje mayor máximo de las elipses en kpc
#define epsilon 0.3                 //Excentricidad de las elipses
#define pi 3.14159265358979323846   //Valor de pi
#define const_k 0.6                 //Constante para el ángulo de rotación de las elipses

//Constantes para el Sol
#define a_Sol 8.6702                //Eje mayor del Sol, cota superior
//#define a_Sol 7.6887              //Eje mayor del Sol, cota inferior

#define epsilon_Sol 0.3             //Excentricidad de la órbita del Sol, común al resto de Estrellas
//#define epsilon_Sol 0.06          //Excentricidad de la órbita del Sol real, para el cálculo del periodo del Sol


///////////////////////////////////////////////////////////////////////////////////////////////////////

//DECLARACION DE FUNCIONES

///////////////////////////////////////////////////////////////////////////////////////////////////////

//Funciones para rellenar matrices y vectores
void rellena_m (double m[]);                                                                        //Función para rellenar el vector m. Todas tendra 1 masa solar
void rellena_inicial_r (double r[][2], double pos_inicial[][3], mt19937_64& generator);             //Función para rellenar la matriz r con las posiciones iniciales de los cuerpos, distribuidos en una espiral logaritmica con ruido aleatorio
void rellena_v(double pos_inicial[][3], double v[][2]);                                             //Función para rellenar la matriz v con las velocidades iniciales, que se calculan a partir de las posiciones iniciales y la masa de Sagitario A*
void rota_elipses(double r[][2], double v[][2], double pos_inicial[][3]);                           //Función para rotar las elipses y con ello las posiciones y velocidades iniciales


//Funciones para trabajar con matrices y vectores
void rellena_cero_matriz (double a[][2]);                                                           //Inicializa a cero una matriz


//Funciones para realizar operaciones matemáticas
void resta_vec (double r[][2], double r_resta[], int v1, int v2);                                   //Función para restar dos vectores
double mod_vec (double v[]);                                                                        //Función para calcular modulo de un vector


//Funciones para realizar el algoritmo
void calculo_a (double r[][2], double a[][2], double m[]);                                          //Función para calcular a
void calculo_r (double r[][2], double v[][2], double a[][2], double h);                             //Funcion para calcular r
void calculo_w (double w[][2], double v[][2], double a[][2], double h);                             //Funcion para calcular w
void calculo_v (double w[][2], double v[][2], double a[][2], double h);                             //Funcion para calcular v


//Funciones para números aleatorios
double r_azar (mt19937_64& generator);                                                              //Función para generar número aleatorio real entre 0 y 1


//Funciones para el Periodo
double calcula_periodo (double r[][2], double r_ant_Sol[], double pos_inicial[][3], double t, bool& periodo_calculado, double& periodo_Sol);




///////////////////////////////////////////////////////////////////////////////////////////////////////

//MAIN

///////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){

    /////////////////////////////////////////////////////////////////////////
    //Definimos variables

    double m[N], r[N][2], r_ant_Sol[2], v[N][2], a[N][2], w[N][2], pos_inicial[N][3];
    bool periodo_calculado = false; //Función para detectar cuando se calcula el periodo
    double periodo_Sol = 0;

    /*
    La matriz pos_inicial se usará para almacenar las posiciones inciales aleatorias de los cuerpos al posiconarlos sobre la elipse. Está formado por el valor del semieje mayor de la elipse y sus coordenadas en polares (r y theta).
    */

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
    rellena_inicial_r(r, pos_inicial, generator);


    //Vector de velocidades iniciales:
    v[0][0] = 0;   //Velocidad inicial de Sagitario A* (cuerpo 0)
    v[0][1] = 0;
    rellena_v(pos_inicial, v);


    //Rotamos las elipses y con ello las posiciones y velocidades iniciales.
    rota_elipses(r, v, pos_inicial);



    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //Abirmos fichero para guardar resultados
    ofstream arch_resultados("resultados_Fase2.txt"); //Abrimos archivo de resultados
    ofstream arch_resultados_simp("resultados_simp_Fase2.txt");   //Abrimos archivo de resultados simplificado
    ofstream arch_cond_iniciales("cond_iniciales_Fase2.txt");   //Abrimos archivo para guardar condiciones iniciales
    ofstream arch_periodo_Sol("periodo_Sol.txt");   //Abrimos archivo para guardar periodo del Sol

    if(!arch_resultados){   //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de resultados" << endl;
    }

    if(!arch_resultados_simp){  //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de resultados" << endl;
    }

    if(!arch_cond_iniciales){  //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de condiciones iniciales" << endl;
    }

    if(!arch_periodo_Sol){  //Verificamos que se haya abierto correctamente
        cerr << "Error al abir el archivo de periodo del Sol" << endl;
    }

    //Encabezado de los archivos
    arch_resultados << "Tiempo" << "\t" << "Cuerpo" << "\t" << "rx" << "\t" << "ry" << "\t" << "vx" << "\t" << "vy" << "\t" << "ax" << "\t" << "ay" << endl;
    arch_resultados_simp << "Tiempo" << "\t" << "Cuerpo" << "\t" << "rx" << "\t" << "ry" << endl;
    arch_cond_iniciales << "Cuerpo" << "\t" << "Eje mayor" << "\t" << "r" << "\t" << "theta" << "\t" << "rx" << "\t" << "ry" << "\t" << "vx" << "\t" << "vy" << endl;
    arch_periodo_Sol << "Periodo del Sol (Myr)" << endl;
    
    for (int i=1; i<N; i++){
        arch_cond_iniciales << i << "\t" << pos_inicial[i][0] << "\t" << pos_inicial[i][1] << "\t" << pos_inicial[i][2] << "\t" << r[i][0] << "\t" << r[i][1] << "\t" << v[i][0] << "\t" << v[i][1] << endl;
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //Algoritmo

    //Inicializamos la aceleración a cero y caluclamos la aceleración por primera vez, el resto de veces será en el bucle
    rellena_cero_matriz(a);
    calculo_a(r,a,m);

    int n=0;
    while (n<iteraciones_max){

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //Guardamos la posición del Sol antes de sobreescribir

        for (int i=0; i<=1; i++){
            r_ant_Sol[i] = r[1][i];
        }


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


        //LLamamos a la función para calcular el periodo del Sol
        periodo_Sol = calcula_periodo(r, r_ant_Sol, pos_inicial, t, periodo_calculado, periodo_Sol);
        if (periodo_calculado){
            arch_periodo_Sol << periodo_Sol << endl;
            periodo_calculado = false; //Reiniciamos la función para que siga calculando periodos cada vez que el Sol complete una órbita
        }

        
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


//Función para rellenar la matriz r con las posiciones iniciales de los cuerpos. Repartimos los cuerpos en elipses concéntricas
void rellena_inicial_r(double r[][2], double pos_inicial[][3], mt19937_64& generator){

    //Para el sol
    pos_inicial[1][0] = a_Sol; //Semieje mayor de la elipse del sol
    pos_inicial[1][2] = 0; //Posición angular del sol en la elipse

    pos_inicial[1][1] = pos_inicial[1][0]*(1-pow(epsilon_Sol,2))/(1 - epsilon_Sol*cos(pos_inicial[1][2]));

    r[1][0] = pos_inicial[1][1]*cos(pos_inicial[1][2]);   //Coordenada x del sol
    r[1][1] = pos_inicial[1][1]*sin(pos_inicial[1][2]);   //Coordenada y del sol

    //Resto de cuerpos
    for (int i=2; i<N; i++){

        //Elegimos al azar el semieje mayor de la elipse y la posición angular que ocupa el cuerpo i en la elipse
        pos_inicial[i][0] = a_min + (a_max - a_min)*r_azar(generator); 
        pos_inicial[i][2] = 2*pi*r_azar(generator); 
        //pos_inicial[i][2] = 0;

        //Calculamos la distancia del cuerpo i al centro de la galaxia y proyectamos en x e y
        pos_inicial[i][1] = pos_inicial[i][0]*(1-pow(epsilon,2))/(1 - epsilon*cos(pos_inicial[i][2]));
        r[i][0] = pos_inicial[i][1]*cos(pos_inicial[i][2]);
        r[i][1] = pos_inicial[i][1]*sin(pos_inicial[i][2]);
    }
}


//Función para rellenar la matriz v con las velocidades iniciales.
void rellena_v(double pos_inicial[][3], double v[][2]){

    double C; //Cte auxiliar

    //Para el Sol
    C = sqrt(G_reescalada*(M_sagitarioA_reescalada + 1)/(pos_inicial[1][0]*(1-pow(epsilon_Sol,2))));
    v[1][0] = -C*sin(pos_inicial[1][2]);
    v[1][1] = C*(cos(pos_inicial[1][2])-epsilon_Sol);

    //Resto de cuerpos
    for (int i=2; i<N; i++){

        C = sqrt(G_reescalada*(M_sagitarioA_reescalada + 1)/(pos_inicial[i][0]*(1-pow(epsilon,2))));

        v[i][0] = -C*sin(pos_inicial[i][2]);
        v[i][1] = C*(cos(pos_inicial[i][2])-epsilon);

    }

    //La distinción entre cuerpos se debe a que usamos una excentricidad distinta para el Sol que para el resto de cuerpos
}


//Funcion para inicializar a cero
void rellena_cero_matriz (double a[][2]){
    for (int i=0; i<N; i++){
        for (int j=0; j<2; j++){
            a[i][j] = 0;
        }
    }
}


//Funcion para rotar las elipses.
void rota_elipses(double r[][2], double v[][2], double pos_inicial[][3]){

    double phi_rotacion, x_rot, y_rot, vx_rot, vy_rot; //Ángulo de rotación

    for (int i=1; i<N; i++){

        if ((i%2)==0) {
            phi_rotacion = const_k * pos_inicial[i][0];
        } 
        else {
            phi_rotacion = const_k * pos_inicial[i][0] + pi;
        }


        //Rotamos la posición
        x_rot = r[i][0]*cos(phi_rotacion) + r[i][1]*sin(phi_rotacion);
        y_rot = -r[i][0]*sin(phi_rotacion) + r[i][1]*cos(phi_rotacion);
        r[i][0] = x_rot;
        r[i][1] = y_rot;

        //Rotamos la velocidad
        vx_rot = v[i][0]*cos(phi_rotacion) + v[i][1]*sin(phi_rotacion);
        vy_rot = -v[i][0]*sin(phi_rotacion) + v[i][1]*cos(phi_rotacion);
        v[i][0] = vx_rot;
        v[i][1] = vy_rot;
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



///////////////////////////////////////////////////////////////////////////////////////////////////////
//Funciones para el periodo del Sol

double calcula_periodo (double r[][2], double r_ant_Sol[], double pos_inicial[][3], double t, bool& periodo_calculado, double& periodo_Sol){

    if (!periodo_calculado){
        
        double y_Sol = r[1][1];
        double y_anterior_Sol = r_ant_Sol[1];

        //Vamos a desrotar las coordenadas del Sol para simplificar el cálculo del periodo.

        double theta = -(pos_inicial[1][0] * const_k + pi); //Ángulo de rotación de las elipses
        double y_Sol_desrot = -r[1][0]*sin(theta) + r[1][1]*cos(theta);
        double y_anterior_Sol_desrot = -r_ant_Sol[0]*sin(theta) + r_ant_Sol[1]*cos(theta);

        if ((y_anterior_Sol_desrot<0) && (y_Sol_desrot>=0)){
            periodo_Sol = t;
            periodo_calculado = true;
            return periodo_Sol;
        }
    }
    return 0; // Devolver 0 si no se ha calculado el periodo
}