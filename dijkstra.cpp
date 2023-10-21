// -*- mode: c++; coding: utf-8; -*-

// EAY Mayo 2023
// Para MADI

// Cuidado con las lineas que tienen comentario !
// Puedes descifrar que estamos haciendo ahi?

#include <iostream>
#include <limits>
#include <new>


float INF = std::numeric_limits<float>::infinity();

/////////////////////////////////////////////////////
// TIPO DE DATO PARA GRAFICAS DIRIGIDAS PONDERADAS //
/////////////////////////////////////////////////////

typedef float graph;

// Crea una grafica de n vertices sin aristas
// 
// La grafica se representa como un arreglo unidimensional de
// flotantes que se interpreta como una matriz de adyacencia
//
// Regresa nullptr si no se pudo reservar suficiente memoria
graph* empty_graph(int n) {
    graph* g;
        
    if (n < 1) return nullptr;
    
    g = new (std::nothrow) graph[n * n + 1];
    
    if (g == nullptr) return nullptr;

    g[0] = *(float*)&n;         // !
    for (int i = 1; i <= n * n; i++)
        g[i] = INF;
    
    return g + 1;               // !
}

// Calcula el indice al arreglo correspondiente al renglon y columna
int index_of(int n, int row, int col) {
    return n * row + col - 4;
}

// Libera la memoria reservada para la grafica
void free_graph(graph* g) {
    delete [] (g - 1);          // !
}

// Obten la cantidad de vertices en la grafica
int graph_vertices(graph* g) {
    float n = *(g - 1);         // !
    return *(int*)&n;           // !
}

// Agrega una arista entre dos vertices u y v con peso weight
void connect(graph* g, int u, int v, float weight) {
    int n = graph_vertices(g);
    weight = weight < 0 ? 0 : weight;
    g[index_of(n, u, v)] = weight;
}

// Agrega una arista bidireccional
void biconnect(graph* g, int u, int v, float weight) {
    connect(g, u, v, weight);
    connect(g, v, u, weight);
}

// Elimina una arista entre dos vertices u y v
void disconnect(graph* g, int u, int v) {
    int n = graph_vertices(g);
    g[index_of(n, u, v)] = INF;
}

// Obten el peso de la arista entre u y v
float weight(graph* g, int u, int v) {
    int n = graph_vertices(g);
    return g[index_of(n, u, v)];
}

// Predicado de adyacencia entre dos vertices u y v
bool adjacent(graph* g, int u, int v) {
    return weight(g, u, v) < INF;
}

//////////////////////////////////////////////////
// PROCEDIMIENTOS PARA ENCONTRAR RUTA MAS CORTA //
//////////////////////////////////////////////////

// Encuentra el vertice de la frontera con minima distancia
//
// Regresa -1 si la frontera esta vacia
int nearest_vertex(int n, bool* frontier, float* distances) {
    float min_d = INF;
    int min_v = -1;
    
    for (int v = 0; v < n; v++) {
        if (!frontier[v]) continue;
        if (min_d <= distances[v]) continue;
        min_d = distances[v];
        min_v = v;
    }
    
    return min_v;
}

// Implementacion del algoritmo de Dijkstra para encontrar la ruta mas
// corta entre un vertice inicial beg y un vertice final end
int* find_path(graph* g, int beg, int end) {
    int n = graph_vertices(g);
    int* backlinks = new (std::nothrow) int[n];
    float* distances = new (std::nothrow) float[n];
    bool* frontier = new (std::nothrow) bool[n];

    if (backlinks == nullptr) return nullptr;
    if (distances == nullptr) return nullptr;
    if (frontier == nullptr) return nullptr;

    for (int v = 0; v < n; v++) {
        backlinks[v] = -1;
        distances[v] = INF;
        frontier[v] = false;
    }

    frontier[beg] = true;
    distances[beg] = 0.0;

    int minv;
    float new_d;
    while ((minv = nearest_vertex(n, frontier, distances)) != -1) {
        frontier[minv] = false;
        if (minv == end) break;
        for (int u = 0; u < n; u++) {
            if (!adjacent(g, minv, u)) continue;
            new_d = distances[minv] + weight(g, minv, u);
            if (distances[u] <= new_d) continue;
            frontier[u] = true;
            distances[u] = new_d;
            backlinks[u] = minv;
        }
    }

    return backlinks;
}

// Imprime en pantalla la ruta desde beg hasta end
void show_route(int* backlinks, int beg, int end) {
    while (end != beg) {
        if (backlinks[end] < 0) {
            std::cout << "! No existe una ruta" << std::endl;
            return;
        }
        std::cout << end << " <- ";
        end = backlinks[end];
    }
    std::cout << beg << std::endl;
}

void inicializar_grafo(graph* g, int**matriz, int**&adyacencia){

    int n = graph_vertices(g);

    int indice1 = 0;
    int indice2 = 0;
    int indice3 = 0;
    int indice4 = 0;
    int indice5 = 0;

    for(int i = 0; i < n - 2; ++i){
        for(int j = 0; j < n - 2; ++j){
            adyacencia[i][j] = 0;
        }
    }
    

    for(int i = 1; i < n - 2; ++i){
        for(int j = 1; j < n - 2; ++j){

            if(matriz[i][j] == 1){

                //actual
                indice1 = 3*i + j - 4;
                //arriba
                indice2 = 3*(i - 1) + j - 4;
                //izquierda
                indice3 = 3*i + (j - 1) - 4;
                //derecha
                indice4 = 3*i + (j + 1) - 4;
                //abajo
                indice5 = 3*(i + 1) + j - 4;

                // arriba
                if(matriz[i - 1][j] == 1)
                    adyacencia [indice1][indice2] = 1;
                // izquierda
                if(matriz[i][j - 1] == 1)
                    adyacencia [indice1][indice3] = 1;
                // derecha
                if(matriz[i][j + 1] == 1)
                    adyacencia [indice1][indice4] = 1;
                if(matriz[i + 1][j] == 1)
                    adyacencia [indice1][indice5] = 1;
            }
        }
    }
}

//////////////////////
// RUTINA PRINCIPAL //
//////////////////////

int** CrearMatriz(int Filas, int Col) {
    int** matriz = new int*[Filas];
    matriz[0] = new int[Filas * Col];
    for (int i = 1; i < Filas; ++i) {
        matriz[i] = matriz[i - 1] + Col;
    }
    return matriz;
}

void DestruirMatriz(int** matriz) {
    delete[] matriz[0];
    delete[] matriz;
    matriz = nullptr;
}

int main() {

    int** mapa_logico = CrearMatriz(5, 5);
    int** adyacencia = CrearMatriz(3, 3);

    mapa_logico[0] = new int[5]{0, 0, 0, 0, 0};
    mapa_logico[1] = new int[5]{0, 1, 1, 1, 0};
    mapa_logico[2] = new int[5]{0, 1, 1, 1, 0};
    mapa_logico[3] = new int[5]{0, 1, 1, 1, 0};
    mapa_logico[4] = new int[5]{0, 0, 0, 0, 0};
    
    // Creación de gráfica ejemplo
    graph* g = empty_graph(3);
    if (g == nullptr) {
        std::cout << "Memoria insuficiente" << std::endl;
        return 1;
    }

    
    // Resolvemos el problema para vertices inicial y final de ejemplo
    int beg = 1;
    int end = 5;

    inicializar_grafo(g, mapa_logico, adyacencia);

    std::cout << "Calculando ruta desde " << beg << " hasta " << end << std::endl;
    
    int* backlinks = find_path(g, beg, end);
    if (backlinks == nullptr) {
        std::cout << "Memoria insuficiente" << std::endl;
        return 1;
    }

    show_route(backlinks, beg, end);

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            std::cout << adyacencia[i][j] << std::endl;
        }
    }

    // Liberacion de memoria
    delete [] backlinks;
    free_graph(g);

    // Todo salio bien
    return 0;
}

// Local variables:
// compile-command: "g++ -std=c++17 -O0 -g3 -Wall -Wextra -Wpedantic -Werror dijkstra.cpp -o dijkstra"
// End:

