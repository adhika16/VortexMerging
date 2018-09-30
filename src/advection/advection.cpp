#include "../../global.hpp"
#include "advection.hpp"

mat advection::advect(const mat& XPOLD, float h, const mat& w, float beta, float gamma){
	mat velocity(XPOLD.n_rows,2);
	vec v1(3);
	vec v2(3);

	for (int i = 0; i < XPOLD.n_rows; i++)
	{
		mat sumU(XPOLD.n_rows,2,fill::zeros);
		for (int j = 0; j < XPOLD.n_rows; j++)
		{
			if( i != j){
				mat dX 		= XPOLD.row(i) - XPOLD.row(j);
				float R 	= norm(dX);
				mat rho 	= dX / h;
				float r 	= norm(rho);
				float q 	= 1 / (2*pi) * (1 - exp(-0.5 * r * r));
				v1(0) 		= dX(0,0);	 
				v1(1)		= dX(0,1);
				v1(2)		= 0;
				v2(0) 		= 0;	 
				v2(1)		= 0;
				v2(2)		= w(j,0) * pow(h,2);
				vec v1v2  	= q / (pow(R,2)) * cross(v1,v2);
				sumU(j,0)	= v1v2(0);
				sumU(j,1)	= v1v2(1);
			}
		}
		float dv_x 		= - beta*XPOLD(i,0);
		float dv_y		= (beta-gamma)*XPOLD(i,1);
		velocity(i,0) 	= - accu(sumU.col(0)) + dv_x;
		velocity(i,1) 	= - accu(sumU.col(1)) + dv_y;
	}

	return velocity;
}
