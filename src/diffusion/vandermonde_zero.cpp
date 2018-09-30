#include "../../global.hpp"
#include "dc_operator.hpp"


mat dc_operator::Vandermonde_zero(mat& Xp, mat& Xq, float eps, int j){
	/*Xp is source points, Xq is collocation points*/
	mat pp;
	mat dx = ( Xp.col(0) - Xq.col(0) ) / eps;
	mat dy = ( Xp.col(1) - Xq.col(1) ) / eps;

	if (j == 0){
		pp = 1;
	}
	else if (j == 1){
		pp = dx;
	}
	else if (j == 2){
		pp = dy;
	}
	else if (j == 3){
		pp = dx*dy;
	}
	else if (j == 4){
		pp = dx*dx;
	}
	else if (j == 5){
		pp = dy*dy;
	}
	else if (j == 6){
		pp = dx*dx*dx;
	}
	else if (j == 7){
		pp = dx*dx*dy;
	}
	else if (j == 8){
		pp = dx*dy*dy;
	}
	else if (j == 9){
		pp = dy*dy*dy;
	}
	else if (j == 10){
		pp = dy*dy*dy*dy;
	}
	else if (j == 11){
		pp = dx*dy*dy*dy;
	}
	else if (j == 12){
		pp = dx*dx*dy*dy;
	}
	else if (j == 13){
		pp = dx*dx*dx*dy;
	}
	else {
		pp = dx*dx*dx*dx;
	}

	return pp;
}
