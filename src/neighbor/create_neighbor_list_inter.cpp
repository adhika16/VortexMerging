#include "../../global.hpp"
#include "neighbor.hpp"


field<mat> neighbor::create_neighbor_list_inter(const mat& dat_fit, const mat& dat_temp, const mat& Rc){
	/*X_temp and Rc_temp are from actual and ghost particle before particle removal-insertion*/
	int sz 			= 	dat_fit.n_rows;
	field<mat> neighbor_list(1,sz);
	for (int i = 0; i < sz; i++)
	{
		mat x_i 			= dat_fit(i, span(0,1));
		mat R_i 			= Rc.row(i);
		neighbor_list(0,i) 	= neighbor_between(x_i, dat_temp, accu(R_i), Rc);
	}
	return neighbor_list;
}
