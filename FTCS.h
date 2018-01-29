

#ifndef FTCS_H_
#define FTCS_H_

/*
 * CrankNicholson.h
 *
 *  Created on: 27 oct. 2017
 *      Author: yoann
 *
 *
 *	This class inherits from Solver and provide an approximated solution of the heat diffusion equation using Dufort & Frankel's explicit scheme
 */

#include "Solver.h"

class FTCS: public Solver {
public:

	//CONSTRUCTOR
	/**
	 * Construcs a solver of the problem using DuFort-Frankel method
	 * @exception invalid_argument ("dx should be positive")
	 * @exception invalid_argument ("dt should be positive")
	 * @exception invalid_argument ("L should be positive")
	 * @exception invalid_argument ("T should be positive")
	 * @exception invalid_argument ("L should be equal or larger than dx")
	 * @exception invalid_argument ("T should be equal or larger than dt")
	 */
	FTCS(double dx/**< double. distance between two space steps */,
			 double dt/**< double. time between two time steps */,
			 double L/**< double. width of the 1D material to consider */,
			 double T/**< double. Total time of the considerated problem */,
			 double D/**< double. Diffusion coefficient of the material */,
			 double Tsur/**< double. The temperature that will be applied on the boundaries of the material */,
			 double Tin/**< double. The initial temperature of the material */,
			 int npes /**<int. The number of processor involved*/,
			 int myRank /**int. The rank of the current processor*/ );
	/**
	 *  Compute the solution and return it. This method is the Dufort-Frankel method applied to the heat diffusion equation problem
	 *  @return Matrix. The computed matrix, can also be accesed through getComputedSolution()
	 *  @see getComputedSolution()
	 */
	Matrix virtual computeSolution();
	virtual ~FTCS();
};

#endif /* FTCS_H_ */
