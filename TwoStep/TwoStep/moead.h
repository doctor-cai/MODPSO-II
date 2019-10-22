#ifndef __MOEAD_H_
#define __MOEAD_H_

#include "scalarfunc.h"
#include "recombination.h"
#include "TSOP.h"
#include <direct.h>
#include<imagehlp.h>

//void save_results(char *archiveName);  /* Write results to file */
//void save_results();
//void save_front();               // save the pareto front into files
//void save_chrom();

class  TMOEAD        
{
public:

	TMOEAD();
	 ~TMOEAD();

	void init_uniformweight();    // initialize the weights for subproblems
	void init_neighbourhood();          // calculate the neighbourhood of each subproblem
	void init_population();             // initialize the population
	void update_weightvector();    // update the weightvector
	void init_velocity();               // initialize the velocity
	void store_pbests(); /* Store personal bests (both variable and fitness values) of particles */
	void update_pbests(int index);
	void update_pbests_child(TIndividual &ind, int index);
	void compute_velocity(int index, TIndividual &child, int gen);
	void update_reference(TIndividual &ind);           // update the approximation of ideal point
//	void update_problem(TIndividual &child, int id, int Tr);   // compare and update the neighboring solutions
	void update_problem(TIndividual &indiv, int id);   // compare and update the neighboring solutions
	void evolution(int gen);                                  // mating restriction, recombination, mutation, update
	void run(int mg, int rn, int ithData, int ithRun);
	void save_results(char savePF[100], char saveChrom[100]);

	
    vector <TSOP>  population;  // current population     
	vector <TIndividual> pbest; // pbest population
	vector <TIndividual> gbest; // gbest population,

	TIndividual *indivpoint;    // reference point
	double **velocity;
	void operator=(const TMOEAD &emo);
};

#endif

