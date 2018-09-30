#include "../../global.hpp"
#include "initialization.hpp"


mat initialization::create_index(const mat& X_particle){
	int sz 		= X_particle.n_rows;
	mat index(sz,1);
	for (int i = 0; i < sz; i++)
	{
		index.row(i) = i;
	}
	return index;
}
