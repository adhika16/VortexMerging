#include "../../global.hpp"
#include "remeshing.hpp"
#include "../../src/diffusion/dc_operator.hpp"


mat remeshing::M4(const mat& dat_uni, const mat& dat_fit_new, const field<mat> neighbor_list_between, float h){
	mat f_interp(dat_uni.n_rows,1);
	vec x_uni 	= dat_uni.col(0);
	vec y_uni 	= dat_uni.col(1);
	vec x_rand 	= dat_fit_new.col(0);
	vec y_rand 	= dat_fit_new.col(1);
	// mat X_uni 	= dat_uni.cols(0,1);
	// mat X_rand 	= dat_fit_new.cols(0,1);
	// vec f_rand 	= dat_fit_new.col(3);
	vec f_rand 	= dat_fit_new.col(2);
	for (int i = 0; i < dat_uni.n_rows; i++)
	{
		mat neighbor = neighbor_list_between(0,i);
		mat W1(neighbor.n_rows,1,fill::zeros);
		mat W2(neighbor.n_rows,1,fill::zeros);
		// mat W(neighbor.n_rows,1,fill::zeros);		
		mat sum_f_interp(neighbor.n_rows,1,fill::zeros);
		for (int j = 0; j < neighbor.n_rows; j++)
		{
			int idx_j		= neighbor(j,0); 	
			float dx 		= abs(x_uni(i) - x_rand(idx_j)) / h;
			float dy 		= abs(y_uni(i) - y_rand(idx_j)) / h;
			float f_tild	= f_rand(idx_j);
			// float r 		= norm(X_uni.row(i) - X_rand.row(idx_j)) / h;
			if (dx >= 1 && dx <= 2 )
			{
				W1(j,0) = 0.5 * pow((2 - dx),2) * (1 - dx);
			}
			if (dx <= 1)
			{
				W1(j,0) = 1 - (5/2)*pow(dx,2) + (3/2)*pow(dx,3);
			}

			if (dy >= 1 && dy <= 2 )
			{
				W2(j,0) = 0.5 * pow((2 - dy),2) * (1 - dy);
			}
			if (dy <= 1)
			{
				W2(j,0) = 1 - (5/2)*pow(dy,2) + (3/2)*pow(dy,3);
			}

			sum_f_interp.row(j) = f_tild * W1(j,0) * W2(j,0);
			// sum_f_interp.row(j) = f_tild * W1(j,0) * W2(j,0); 
		}

		f_interp.row(i) = accu(sum_f_interp);
	}

	return f_interp;
}
