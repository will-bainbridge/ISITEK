#ifndef QUADRATURE_H
#define QUADRATURE_H

double quadrature_gauss_location(int order, int index);
double quadrature_gauss_weight(int order, int index);
double quadrature_hammer_location(int order, int dimension, int index);
double quadrature_hammer_weight(int order, int index);

#endif
