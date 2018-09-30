#ifndef _remeshing_hpp_
#define _remeshing_hpp_

class remeshing
{
public:
	arma::mat dc_pse_interp(const arma::mat& dat_fit_new, const arma::mat& dat_temp, const arma::mat& f_sour, const arma::mat& f_coll, int l, float eps, const arma::mat& Rcp, const arma::field<arma::mat>& neighbor_list /*, int choice*/);
	arma::mat M4(const arma::mat& dat_uni, const arma::mat& dat_fit_new, const arma::field<arma::mat> neighbor_list_between, float h);
};

#endif