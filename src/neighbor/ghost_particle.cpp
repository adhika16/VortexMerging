#include "../../global.hpp"
#include "neighbor.hpp"

mat neighbor::ghost_particle(const mat& XP, const mat& index, float maxdom, float mindom, int choice){
	// dat2ghost = join_rows(XP, index)
	/*
	_________________________
	|		|		|		|		
	|	A	|	B	|	C	|
	|_______|_______|_______|
	|		|		|		|
	|	D	| domain|	E	|
	|_______|_______|_______|
	|		|		|		|
	|	F	|	G	|	H	|
	|_______|_______|_______|
	*/
	double delta 	= 1e-4;
	mat dat2ghost	= join_rows(XP, index);
	int sz 		= dat2ghost.n_rows;
	mat domain	= dat2ghost;
	mat datA	= dat2ghost; 	
	mat datB	= dat2ghost; 
	mat datC	= dat2ghost; 
	mat datD	= dat2ghost; 
	mat datE	= dat2ghost; 
	mat datF	= dat2ghost; 
	mat datG	= dat2ghost; 
	mat datH	= dat2ghost; 

	mat Lxy 	= ones(1,2) * maxdom;
	mat Lyx(1,2); 	Lyx(0,0) 	= maxdom; 	Lyx(0,1) 	= -maxdom;
	mat Lx(1,2); 	Lx(0,0) 	= maxdom; 	Lx(0,1) 	= 0;
	mat Ly(1,2); 	Ly(0,0) 	= 0; 		Ly(0,1) 	= maxdom;	

	for (int i = 0; i < sz; i++)
	{
		datA(i, span(0,1)) = datA(i, span(0,1)) + (-Lyx) * 2 - delta;
		datB(i, span(0,1)) = datB(i, span(0,1)) + (Ly) * 2 + delta;
		datC(i, span(0,1)) = datC(i, span(0,1)) + (Lxy) * 2 + delta;
		datD(i, span(0,1)) = datD(i, span(0,1)) + (-Lx) * 2 - delta;
		datE(i, span(0,1)) = datE(i, span(0,1)) + (Lx) * 2 + delta;
		datF(i, span(0,1)) = datF(i, span(0,1)) + (-Lxy) * 2 - delta;
		datG(i, span(0,1)) = datG(i, span(0,1)) + (-Ly) * 2 - delta;
		datH(i, span(0,1)) = datH(i, span(0,1)) + (Lyx) * 2 + delta;
	}

	mat matghost 	= join_cols(domain, datA);
	matghost 		= join_cols(matghost, datB);
	matghost 		= join_cols(matghost, datC);
	matghost 		= join_cols(matghost, datD);
	matghost 		= join_cols(matghost, datE);
	matghost 		= join_cols(matghost, datF);
	matghost 		= join_cols(matghost, datG);
	matghost 		= join_cols(matghost, datH);

	sz 				= matghost.n_rows;
	double offs 	= 1.0;
	mat dat_temp 	= zeros(1,3);
	mat dat_fit 	= zeros(1,3);
	for (int i = 0; i < sz; i++)
	{
		float x 	= matghost(i,0);
		float y 	= matghost(i,1);
		if ((x) >= mindom && (x) <= maxdom && (y) >= mindom && (y) <= maxdom ){
			dat_fit.insert_rows(0, matghost.row(i));
		}
		if ((x) >= offs*mindom && (x) <= offs*maxdom && (y) >= offs*mindom && (y) <= offs*maxdom){
			dat_temp.insert_rows(0, matghost.row(i));
		}
	}
	int n_fit 		= dat_fit.n_rows;
	dat_fit.shed_row(n_fit-1);
	int n_temp 		= dat_temp.n_rows;
	dat_temp.shed_row(n_temp-1);	

	if (choice == 1){
		return dat_fit;
	}
	if (choice == 2){
		return dat_temp;
	}
	if (choice == 3){
		return matghost;
	}	
}
