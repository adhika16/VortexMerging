#include "../../global.hpp"
#include "neighbor.hpp"

mat neighbor::neighbor_temp(const mat& x_i, const mat& dat_temp, float Rc_i, const mat& Rcq){
	float R, Rcpq, Rcpq0;
	int nIter 		= 	0;
	int nidx 		= 	0;
	mat Rcp 		= 	ones(1,1) * Rc_i;
	int sz 			= 	dat_temp.n_rows;
	mat xq 			= 	dat_temp.cols( 0,1 );
	vec temp_index	= 	dat_temp.col( 2 );
	mat idx 		= 	zeros( sz,1 );
	for (int i = 0; i < sz; i++) 
	{
		R 			= 	norm ( x_i - xq.row(i) ) ;
		mat Rcq_i 	= 	Rcq.row( temp_index(i) ) ;
		mat RcpRcq 	= 	join_rows( Rcq_i, Rcp ) ;
		Rcpq 		= 	RcpRcq.min();
		if ( R <= Rcpq*(1 + 0.1) && R > 0 ) {
			idx.row(nidx) 	= ( i );
			nidx 			+= 	1;
		}
		else{
			continue;
		}
	}
	int strt 	= 	0;
	int fnsh 	= 	nidx-1;
	mat idx_ 	= 	idx.rows( strt, fnsh );
	return idx_;
}
