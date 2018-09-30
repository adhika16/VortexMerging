#ifndef _dc_operator_hpp_
#define _dc_operator_hpp_


class dc_operator
{
	arma::mat Vandermonde(arma::mat &Xp, arma::mat &Xq, float eps, int j);
	arma::mat Vandermonde_zero(arma::mat &Xp, arma::mat &Xq, float eps, int j);
public:
	arma::mat Q22(const arma::mat& dat_fit, const arma::mat& dat_temp, const arma::mat& matf, const arma::mat& matf_old, 
		float beta, int l, double eps, const arma::mat& Rcp, const arma::field<arma::mat>& neighbor_list);

	friend class remeshing;
};


#endif