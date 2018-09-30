#include "../../global.hpp"
#include "initialization.hpp"


mat initialization::init_uniform(int n, float h, float mindom){
	mat X_uniform (n * n,2);
	int it = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			X_uniform(it,0) 	= ( (i+1) - 0.5 ) * h + mindom;
			X_uniform(it,1) 	= ( (j+1) - 0.5 ) * h + mindom;
			it 					+= 1;
		}
	}

	return X_uniform;
}
