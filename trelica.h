#include "aux.h"

#pragma once

void item_c();


struct no{
    double massa;
    double x, y;
};
typedef struct no no;

struct barra {
    no* no_1;
    no* no_2;
    double angulo;
    double comprimento;
};
typedef struct barra barra;

struct trelica{
    int n_barras;
    int n_nos;
    int n_nos_nao_fixos;
    int ro, E;
    double A;

    barra** barras;
    no** nos;
};

typedef struct trelica trelica;
