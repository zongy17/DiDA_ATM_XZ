//
// Created by lenovo on 2020/12/3.
//

#ifndef XZ_V11_PARAMS_H
#define XZ_V11_PARAMS_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <armadillo>
using namespace std;
using namespace arma;

typedef double type_f;

#define sound_velocity 340
#define ph0 1e5
#define p0 1e5
#define Rd 287.04
#define gamma 1.4
#define gravity 9.80616
#define pi 3.1415926

void Thomas_solve(Col<type_f> & A, Col<type_f> & B, Col<type_f> & C, Col<type_f> & D, Col<type_f> & X);

#endif //XZ_V11_PARAMS_H
