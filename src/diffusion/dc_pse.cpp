#include "../../global.hpp"
#include "dc_operator.hpp"

mat dc_operator::Q22(const mat& dat_fit, const mat& dat_temp, const mat& matf, const mat& matf_old, float beta, int l, double eps, 
	const mat& Rcp, const field<mat>& neighbor_list){
	/*	matf is matrix containing a continuous function to be analyze */
	// Xsour : XPOLD or dat_temp
	// Xcoll : XPNEW or dat_fit
	int sz 			= dat_fit.n_rows;
	double 	Rc; // eps;
	mat 	Q22_dc(sz,1);
	mat Xcoll 		= dat_fit.cols(0,1);
	mat Xsour 		= dat_temp.cols(0,1);
	vec temp_index	= dat_temp.col(2);
	vec fit_index	= dat_fit.col(2);

	/*looping started ...*/
	for (int K1 = 0; K1 < sz; K1++)
	{
		Rc 			= (Rcp ( fit_index(K1),0 ) );
		// eps 		= Rc;
		mat idx 	= neighbor_list(0,K1);
		int nidx 	= idx.n_rows;
	// arrange matrix E(xp) 
		mat E 		= zeros(nidx,nidx);
		for (int i = 0; i < nidx; i++)
		{
			mat xp 		= Xcoll.row(K1);
			int idx_i 	= idx(i,0);
			mat xq 		= Xsour.row(idx_i);
			double R 	= norm(xp-xq);
			E(i,i) 		= ( exp( -pow(R/eps,2) * 0.5 ) );
		}
	// arrange matrix V(xp)
		mat V = zeros(nidx,l);
		for (int i = 0; i < nidx; i++)
		{
			for (int j = 0; j < l; j++)
			{
				mat xp 		= Xcoll.row(K1);
				int idx_i 	= idx(i,0);
				mat xq 		= Xsour.row(idx_i);
				double R 	= norm(xp-xq);
				V(i,j) 		= accu(Vandermonde(xp, xq, eps, j)) ; 
			}
		}
	// arrange matrix B(xp)
		mat B 		= E.t() * V;
	// arrange matrix A(xp)
		mat A 		= B.t() * B;
	// arrange matrix b(xp)
		mat b 		= zeros(l,1);
		b.row(2) 	= 2;
		b.row(3) 	= 2;
	// calculate a(xp) 
		// use Choleski decomposition, 'h' must be less than 1 
		mat L 	= chol(A, "lower");
		mat Li 	= solve(L,b);
		mat a 	= solve(L.t(),Li);
		// mat a 	= solve(A,b);
		mat sumQ22(nidx,1);

		for (int K2 = 0; K2 < nidx; K2++)
		{
			mat xp_k2 	= Xcoll.row(K1);
			int idx_k2 	= norm(idx.row(K2));
			mat xq_k2 	= Xsour.row(idx_k2);
			mat dx_k2 	= xp_k2 - xq_k2;
			double R 	= norm(dx_k2);
		// arrange K_epsilon(xp)
			double K_e 	= 0;
			mat K_e_i(l,1);
			for (int i = 0; i < l; i++)
			{
				K_e_i.row(i) 	= a.row(i) * Vandermonde(xp_k2, xq_k2, eps, i);
			}
			K_e 				= accu(K_e_i);
		// arrange phi_e(xp)
			double phi_e 		= exp( -pow( R/eps , 2) );
		// arrange eta_epsilon^beta(xp)
			double eta_e 		= phi_e * K_e;
		// calculate Q(2,2)_f(xp)
			mat dw 				= matf_old.row(temp_index(idx_k2)) - matf.row(fit_index(K1));
			sumQ22.row(K2) 		= dw * eta_e; 
		}
		Q22_dc.row(K1) 			= accu(sumQ22) / pow(eps,2);
	}
	return Q22_dc;
}
