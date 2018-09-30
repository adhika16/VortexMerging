#include "../../global.hpp"
#include "remeshing.hpp"
#include "../../src/diffusion/dc_operator.hpp"

mat remeshing::dc_pse_interp(const mat& dat_fit_new, const mat& dat_temp, const mat& f_sour, const mat& f_coll, int l, float eps, const mat& Rcp, const field<mat>& neighbor_list /*, int choice*/){
	// calling friend function
	dc_operator DC;

	/*matG 			: old particle intensities, 
	Rcp 			: new particle cutoff radius,
	neihgbor_list 	: neihgbor list between XPNEW and XPOLD*/
	int sz 			= 	dat_fit_new.n_rows;
	vec fit_index	= 	dat_fit_new.col(2);
	vec temp_index	= 	dat_temp.col(2);
	l 				= 	l + 1;
	mat XPNEW 		= 	dat_fit_new.cols(0,1);
	mat XPOLD 		= 	dat_temp.cols(0,1);
	mat	f_interp(sz,1);
	// looping started ...
	for (int K1 = 0; K1 < sz; K1++)
	{
		mat idx 		= neighbor_list(0,K1);
		int nidx 		= idx.n_rows;
	// arrange matrix E(xp) 
		mat E 			= zeros(nidx,nidx);
		for (int i = 0; i < nidx; i++)
		{
			mat xp 		= XPNEW.row(K1);
			int idx_i 	= idx(i,0);
			mat xq 		= XPOLD.row(idx_i);
			float R 	= norm(xp-xq);
			E(i,i) 		= pow( (pow((2*R/eps), 14) + 1) , -0.5) ;
		}
	// arrange matrix V(xp)
		mat V 			= zeros(nidx,l);
		for (int i = 0; i < nidx; i++)
		{
			for (int j = 0; j < l; j++)
			{
				mat xp 		= XPNEW.row(K1);
				int idx_i 	= idx(i,0);
				mat xq 		= XPOLD.row(idx_i);
				float R 	= norm(xp-xq);
				V(i,j) 		= accu ( DC.Vandermonde_zero ( xp, xq, eps, j ) ) ; 
			}
		}
	// arrange matrix B(xp)
		mat B 	= 	E.t() * V;
	// arrange matrix A(xp)
		mat A 	= 	B.t() * B;
	// arrange matrix b(xp)
		mat b 	= 	zeros(l,1);
		b(0,0) 	= 	1;	
	// calculate a(xp) 
		// use Choleski decomposition, 'h' must be less than 1 
		// mat L 	= 	chol(A, "lower");
		// mat L2	= 	solve(L,b);
		// mat a 	= 	solve(L.t(), L2);
		mat a  	= 	solve(A,b);
		mat sum_f_interp	= 	zeros(nidx,1);

		for (int K2 = 0; K2 < nidx; K2++)
		{
			mat xp_k2 	= XPNEW.row(K1);
			int idx_k2 	= idx(K2,0);
			mat xq_k2 	= XPOLD.row(idx_k2);
			mat dx_k2 	= xp_k2 - xq_k2;
			double R 	= norm(dx_k2);
		// arrange K_epsilon(xp)
			double K_e 	= 	0;
			mat K_e_i	= 	zeros(l,1);
			for (int i = 0; i < l; i++)
			{
				K_e_i.row(i) 	= a.row(i) * DC.Vandermonde_zero(xp_k2, xq_k2, eps, i);
			}
			K_e 		= 	accu(K_e_i);
		// arrange phi_e(xp)
			double phi_e 			= pow( (pow((2*R/eps), 14) + 1) , -1 );
		// arrange eta_epsilon^beta(xp)
			double eta_e 			= phi_e * K_e;
		// calculate sum_f_interp(xp)
			sum_f_interp.row(K2) 	= f_sour.row( (idx_k2) ) * eta_e; 
			// sum_Q_interp.row(K2) 	= ( matG((idx_k2),0) - mat_f_new((K1),0) ) * eta_e;
		}
		f_interp(K1,0) 				= accu(sum_f_interp);
		// Q_interp(K1,0) 				= accu(sum_Q_interp) / pow(eps,2);
	}

	return f_interp;
}