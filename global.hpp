#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <string>

using namespace std;
using namespace arma;

const float pi = 3.141592654;

// domain parameter
const float maxdom 	= 7.5;
const float mindom 	= -7.5;
const float dev 		= 0.25; // [0.2, 0.25, 0.5]
const float h 		= 0.2;
const float c 		= 1.0; // [0.5, 0.9, 1.4]

// DC PSE operators parameter
const float beta 		= 	2;
const int l_2_0_2		= 	9;
const int l_2_0_4 	= 	20;
const int l_1_0_2		= 	6;
const int l_1_0_4		= 	15;
const int nIter 		= 	0;
// time step parameter
const float ti 		= 	0.0;
const float tf 		= 	1.0;
const float dt 		= 	0.001;
const int N 			= 	32;
const float gammas 	= 	4;
const float betha 	= 	2;
const float nu		= 	1;
/*float dt		= 	1 / (N * Re);*/
// input for Gauss pulse
const float Re	= 10;


#endif
