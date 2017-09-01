#pragma once

#include <stdlib.h>
#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include "armadillo"

#ifdef __APPLE__
    #define _TERMINAL "wxt"
#else
    #define _TERMINAL "wxt"
#endif

using namespace std;
using namespace arma;

FILE *gpc_init_image();
int gpc_plot_image(FILE *pipe, vector<vec> vectors_A, vec w, vec alpha_vec, vec chi2_vec, double lalpha_min, int indexToPlot, int bosonOffset);
void gpc_close(FILE *pipe);

