#include "../../global.hpp"
#include "functions.hpp"

mat functions::f_merger_vortex(const mat& XPOLD, const mat& X0, float Re){
	// test function 
	float q 	= 2.56085;
	float A 	= 20;
	float R0 	= 0.8;
	float Gamma = 2 * pi * Re;
	int sz 		= XPOLD.n_rows;
	vec x 		= XPOLD.col(0);
	vec y 		= XPOLD.col(1);

	float x0 	= X0(0,0);
	float y0	= X0(0,1);
	mat f(sz,1,fill::zeros);
	for (int i = 0; i < sz; i++) {
		float d1 	= - pow((x(i) - x0), 2) - pow(y(i), 2);
		float d2 	= - pow((x(i) + x0), 2) - pow(y(i), 2);
		f(i,0)		= - Gamma / (2 * pi) * ( exp(d1) + exp(d2) );
	}
	return f;
}