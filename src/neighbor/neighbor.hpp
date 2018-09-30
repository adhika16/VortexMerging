#ifndef _neighbor_hpp_
#define _neighbor_hpp_

class neighbor
{
	arma::mat neighbor_temp (const arma::mat& x_i, const arma::mat& dat_temp, float Rc_i, const arma::mat& Rcq);
	arma::mat neighbor_between(const arma::mat& x_i, const arma::mat& dat_temp, float Rc_i, const arma::mat& Rcq);
public:
	arma::mat ghost_particle_unsym_v2(const arma::mat& XP, const arma::mat& index, const arma::mat& f, const arma::mat& u, float maxdom, float mindom, int choice);
	arma::mat ghost_particle_unsym(const arma::mat& XP, const arma::mat& index, float maxdom, float mindom, int choice);
	arma::mat ghost_particle (const arma::mat& XP, const arma::mat& index, float maxdom, float mindom, int choice);
	arma::field<arma::mat> create_neighbor_list (const arma::mat& dat_fit, const arma::mat& dat_temp, const arma::mat& Rc);
	arma::field<arma::mat> create_neighbor_list_inter(const arma::mat& dat_fit, const arma::mat& dat_temp, const arma::mat& Rc);
};


#endif