#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <algorithm> /* implementa��es de algoritmos de ordena��o, soma acumulada, revers�o de vetores e m�ximo e m�nimos de containers */
#include <cmath> /* implementa��es de algoritmos matem�ticos e estat�sticos */
#include <fstream> /* manipula��o de arquivos */
#include <iomanip> /* manipula��o de sa�da e entrada de dados */
#include <iostream> /* implementa��es b�sicas da linguagem */
#include <utility> /* implementa��es de alguns containers e opera��es de swap */
#include <vector> /* implementa��o do container vector e suas opera��es */
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue> // Para a busca em largura
#include <limits> // Para infinito
#include <set> // Para o algoritmo de Kruskal
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <map>
#include <chrono>
#define DBL_MAX __DBL_MAX__
#include <stdexcept>
#include <chrono>
#include <memory> 
#include <cfloat>
#include <numeric>
#include <list>
#include <random>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sys/timeb.h>
//#include <sys/resource.h>
#include <unistd.h>
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL	/* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector */
static int mti = N + 1;		/* mti==N+1 means mt[N] is not initialized */

#endif  //DEFINES_HPP
