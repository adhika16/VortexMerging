#include "../../global.hpp"
#include "initialization.hpp"


mat initialization::init_random(float dev, float mindom, float h, int sz){
	mat x_random(sz * sz,2);
	int it = 0;
	for (int i = 0; i < sz; i++)
	{
		for (int j = 0; j < sz; j++)
		{
			x_random(it,0) = ( (i+1) - 0.5 + -dev + 2 * dev * randu() ) * h + mindom;
			x_random(it,1) = ( (j+1) - 0.5 + -dev + 2 * dev * randu() ) * h + mindom; 
			it += 1;
		}
	}

	return x_random;
}
