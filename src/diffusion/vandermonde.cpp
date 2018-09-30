#include "../../global.hpp"
#include "dc_operator.hpp"

mat dc_operator::Vandermonde(mat& Xp, mat& Xq, float eps, int j){
	mat pp;
	mat dx = ( Xp.col(0) - Xq.col(0) ) / eps;
	mat dy = ( Xp.col(1) - Xq.col(1) ) / eps;

	if (j == 0){
		pp = dx;
	}
	else if (j == 1){
		pp = dy;
	}
	else if (j == 2){
		pp = dx*dx;
	}
	else if (j == 3){
		pp = dy*dy;
	}
	else if (j == 4){
		pp = dx*dy;
	}	
	else if (j == 5){
		pp = (dx*dx)*dx;
	}
	else if (j == 6){
		pp = (dy*dy)*dy;
	}
	else if (j == 7){
		pp = (dx*dx)*dy;
	}
	else if (j == 8){
		pp = dx*(dy*dy);
	}
	else if (j == 9){
		pp = dy*dy*dy*dy;
	}
	else if (j == 10){
		pp = dx*dy*dy*dy;
	}
	else if (j == 11){
		pp = dx*dx*dy*dy;
	}
	else if (j == 12){
		pp = dx*dx*dx*dy;
	}
	else if (j == 13){
		pp = dx*dx*dx*dx;
	}
	else if (j == 14){
		pp = dy*dy*dy*dy*dy;
	}
	else if (j == 15){
		pp = dx*dy*dy*dy*dy;
	}
	else if (j == 16){
		pp = dx*dx*dy*dy*dy;
	}
	else if (j == 17){
		pp = dx*dx*dx*dy*dy;
	}
	else if (j == 18){
		pp = dx*dx*dx*dx*dy;
	}
	else {
		pp = dx*dx*dx*dx*dx;
	}

	return pp;
}
