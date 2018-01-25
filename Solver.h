/*
 * Solver.h
 *
 *  Created on: 25 oct. 2017
 *      Author: yoann
 *
 *   This class is the base class for classes that whish to implements a solution of the heat diffusion equation problem.
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "matrix.h"

class Solver{

protected:
	Matrix computedSolution;
	double dx;
	double dt;
	double L;
	double T;
	double D;
	double Tsur;
	double Tin;
	int npes;
	int myRank;



public:

	//CONSTRUCTOR
	/**
	 * Construcs a solver of the problem, can not be instanciated as an object since it is a virtual base class
	 * @exception invalid_argument ("dx should be positive")
	 * @exception invalid_argument ("dt should be positive")
	 * @exception invalid_argument ("L should be positive")
	 * @exception invalid_argument ("T should be positive")
	 * @exception invalid_argument ("L should be equal or larger than dx")
	 * @exception invalid_argument ("T should be equal or larger than dt")
	 */
	Solver(double dx/**< double. distance between two space steps */,
			 double dt/**< double. time between two time steps */,
			 double L/**< double. width of the 1D material to consider */,
			 double T/**< double. Total time of the considerated problem */,
			 double D/**< double. Diffusion coefficient of the material */,
			 double Tsur/**< double. The temperature that will be applied on the boundaries of the material */,
			 double Tin/**< double. The initial temperature of the material */,
			 int npes /**<int. The number of processor involved*/,
			 int myRank /**int. The rank of the current processor*/);
	//Compute Methods
	/**
	 * Return the computed matrix, computeSolution() has to be called first to get the solution matrix.
	 * @return Matrix, the computed matrix
	 */
	Matrix getComputedSolution() ;

	//ACCESOR METHODS
	/**
	 * get the time step
	 * * @return double. the time step considered
	 */
	double getDT();
	/**
	 * get the space step
	 * * @return double. the space step considered
	 */
	double getDX();
	/**
	 * get the width of the material
	 * * @return double. the width considered
	 */
	double getL();
	/**
	 * get the overall time
	 * * @return double. the overall time considered
	 */
	double getT();
	/**
	 * get the diffusion coefficient
	 * * @return double. the diffusion coefficient considered
	 */
	double getD();
	/**
	 * get the temperature at the surface
	 * * @return double. the temperature applied on the surface
	 */
	double getTsur();
	/**
	 * get the initial temperature
	 * * @return double. the initial temperature considered
	 */
	double getTin();


	//COMPUTATION METHOD

	/**
	 *  Compute the solution and return it. This method must be implemented in the child class if you want it not to be virtual
	 *  @return Matrix. The computed matrix, can also be accesed through getComputedSolution()
	 *  @see getComputedSolution()
	 */

	virtual Matrix computeSolution() = 0;


	//SCREEN/FILE AND OUTPUT
	/**
	 * redifinition of the << operator to the screen, displays the time every 0.1seconde from 0 to 0.5 seconde.
	 * @return the generated stream
	 */
	friend std::ostream& operator<<(std::ostream& os, Solver& m );
	/**
	 * redifinition of the << operator to a file, displays the time every 0.1seconde from 0 to 0.5 seconde. Especially mafe for GNUPlot usage.
	 * @return the generated stream
	 */
	friend std::ofstream& operator<< (std::ofstream& ifs, Solver& m );

	//DESTROYER
	/**
	 * Destroys the object
	 */
	virtual ~Solver(){}
};
#endif /* Solver_h */
