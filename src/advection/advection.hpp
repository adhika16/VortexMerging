#ifndef _advection_hpp_
#define _advection_hpp_

class advection
{
public:
	arma::mat advect(const arma::mat& XPOLD, float h, const arma::mat& w, float beta, float gamma);
};


#endif