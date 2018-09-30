#include "../../global.hpp"
#include "neighbor.hpp"

field<mat> neighbor::create_neighbor_list(const mat& dat_fit, const mat& dat_temp, const mat& Rc){ 
	int sz 			= dat_fit.n_rows;
	field<mat> neighbor_list(1,sz);
	vec fit_index 	= dat_fit.col(2);
	vec temp_index 	= dat_temp.col(2);
	for (int i = 0; i < sz; i++)
	{
		mat x_i 			= dat_fit(i, span(0,1));
		mat R_i 			= Rc.row(fit_index(i));
		neighbor_list(0,i) 	= neighbor_temp(x_i, dat_temp, accu(R_i), Rc);
	}
	return neighbor_list;
}
