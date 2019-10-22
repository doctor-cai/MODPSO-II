#ifndef __TSOP_H_
#define __TSOP_H_

#include "individual.h"

class TSOP 
{
public:
	TSOP();
	virtual ~TSOP();

	TIndividual     indiv;
	vector <double> namda;		// weightvector ,used to the population(÷÷»∫)
	vector <int>    table;     // the vector for the indexes of neighboring subproblems
	vector <int>    tabler;    // Tr replacement neighborhood
	vector <int>    array;
	int n;

    void  operator=(const TSOP&sub2);
};

#endif
