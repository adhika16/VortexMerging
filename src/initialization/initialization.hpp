#ifndef _initialization_hpp_
#define _initialization_hpp_

class initialization
{
public:
	arma::mat init_uniform(int n, float h, float mindom);
	arma::mat init_random(float dev, float mindom, float h, int sz);
	arma::mat create_index(const arma::mat &XP);
};

#endif