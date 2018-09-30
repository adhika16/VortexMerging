#include "global.hpp"
#include "src/initialization/initialization.hpp"
#include "src/initialization/functions.hpp"
#include "src/neighbor/neighbor.hpp"
#include "src/diffusion/dc_operator.hpp"
#include "src/advection/advection.hpp"
#include "src/remeshing/remeshing.hpp"

// =================================================================================== 
int main(int argc, char const *argv[])
{
    // creating instance ============================================================= 
    initialization    init; 
    functions         test_function;
    neighbor          neighboring;
    dc_operator       DC;
    advection 		  adv;
    remeshing 		  remesh;
	// dictionary ===================================================================== 
	int n,n2,nIter;
	float sig,epsilon;
	float t = ti;	
	// domain parameter ===============================================================
	mat X0(1,2); /*center of pulse location*/ 
	X0(0,0)	= 2.0; /*x*/  
	X0(0,1)	= 0.0; /*y*/
	n 		= std::ceil((maxdom-mindom)/h); 
	epsilon = h/c;
	n2		= n*n;
	// initializing domain ============================================================
	mat Rc 				= 	ones(n*n,1) * 3.5 * epsilon; /*radius of influence*/
	mat X_UNI 		 	= 	init.init_uniform(n,h,mindom);
	mat X_RAND		 	= 	init.init_random(dev,mindom,h,n);
	mat XPOLD 		 	= 	X_UNI;
	mat index 		 	= 	init.create_index(XPOLD);
	mat fa			 	= 	test_function.f_merger_vortex(XPOLD, X0, Re); // actual function field
	mat u 			 	= 	adv.advect(XPOLD, h, fa, betha, gammas);
	mat dat_fit			= 	neighboring.ghost_particle_unsym_v2(XPOLD, index, fa, u, maxdom, mindom, 1);
	mat dat_temp		= 	neighboring.ghost_particle_unsym_v2(XPOLD, index, fa, u, maxdom, mindom, 2);
	mat dat_uni 	 	= 	dat_fit;
	mat dat_temp_old 	= 	dat_temp;
	// create neighbor list
	field<mat> neighbor_list = neighboring.create_neighbor_list(dat_fit, dat_temp, Rc); // output : particle indices
	field<mat> neighbor_list_between(1,dat_fit.n_rows);

// ITERATIVE ALGORITHM ================================================================
	cout << " >>> Iteration Begin! ";
	do {
		// convection =================================================================
		mat up 		= adv.advect(dat_fit.cols(0,1), h, dat_fit.col(3), betha, gammas);
		// diffusion ==================================================================
		mat Q22 	= DC.Q22(dat_fit, dat_temp, dat_fit.col(3), dat_temp.col(3), beta, 
						l_2_0_2, epsilon, Rc, neighbor_list);

		// Runge-Kutta 2nd order, time stepping =======================================
		mat xp1		= dat_fit.cols(0,1) + dt * up;
		mat fp1		= dat_fit.col(3) + dt * ( gammas * dat_fit.col(3) + nu * Q22 );
		dat_fit.cols(0,1)	= xp1;
		dat_fit.col(3)		= fp1;
		// create new neighboring
		mat dat_fit_new		= neighboring.ghost_particle_unsym_v2(dat_fit.cols(0,1), dat_fit.col(2), 
								dat_fit.col(3), u, maxdom, mindom, 1);
		mat dat_temp_new	= neighboring.ghost_particle_unsym_v2(dat_fit.cols(0,1), dat_fit.col(2), 
								dat_fit.col(3), u, maxdom, mindom, 2);
		neighbor_list 		= neighboring.create_neighbor_list(dat_fit_new, dat_temp_new, Rc); // output : particle indices

		// remeshing ===================================================================
		if (nIter%5 == 0)
		{
			neighbor_list_between 	= neighboring.create_neighbor_list(dat_uni, dat_fit_new, Rc); // output : particle indices
			mat f_uni 				= remesh.dc_pse_interp(dat_uni, dat_fit_new, dat_fit_new.col(3), 
										dat_fit_new.col(3), l_2_0_2, epsilon, Rc, neighbor_list_between);
			dat_uni.col(3) 			= f_uni;
			dat_temp_old 			= neighboring.ghost_particle_unsym_v2(dat_uni.cols(0,1), dat_uni.col(2), 
										dat_uni.col(3), dat_uni.cols(4,5), maxdom, mindom, 2);
			neighbor_list 			= neighboring.create_neighbor_list(dat_uni, dat_temp_old, Rc); 
			dat_fit 				= dat_uni;
			dat_temp 				= dat_temp_old;
			cout << "Y";			
		}
		else
		{
			dat_fit 	= dat_fit_new;
			dat_temp 	= dat_temp_new;
		}
 		// looping =====================================================================
 		cout << "\nt = " << t;
 		t 		+= dt;
 		nIter 	+= 1;
	} while(t<=tf);
	cout << "\n End of iterations <<<\n ...Program End\n";


	return 0;
}

// Notes : ============================================================
/* there are a lot of mess inside this code. I decided to keep this as 
it is. In addition, I will re-create this code after I learn how did a 
particles is added during a simulation.*/ 
// ====================================================================