#include "../global.hpp"
#include "../src/initialization/functions.hpp"
#include "../src/initialization/initialization.hpp"
#include "../src/neighbor/neighbor.hpp"
#include "../src/advection/advection.hpp"
#include "main_init.h"

void main_init::init_part(arma::mat &Rc, arma::mat &X_UNI,
					arma::mat &X_RAND, arma::mat &XPOLD,
					arma::mat &index, arma::mat &fa,
					arma::mat &u, arma::mat &dat_fit,
					arma::mat &dat_temp, arma::mat &dat_uni,
					arma::mat &dat_temp_old, 
					arma::field<arma::mat> &neighbor_list,
					arma::field<arma::mat> &neighbor_list_between)
{
	// create instance 
    initialization    init; 
    functions         test_function;
    neighbor          neighboring;
    advection 		  adv;



	Rc = Rc_temp;
	X_UNI = X_UNI_temp;
	X_RAND = X_RAND_temp;
	XPOLD = XPOLD_temp;
	index = index_temp;
	fa = fa_temp;
	u = u_temp;
	dat_fit = dat_fit_temp;
	dat_temp = dat_temp_temp;
	dat_uni = dat_uni_temp;
	dat_temp_old = dat_temp_old_temp;
	neighbor_list = neighbor_list_temp;
	neighbor_list_between = neighbor_list_between_temp;
}