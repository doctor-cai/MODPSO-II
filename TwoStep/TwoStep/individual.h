#ifndef __TINDIVIDUAL_H_
#define __TINDIVIDUAL_H_

//#include "global.h"

#include "objectives.h"

class TIndividual
{
public:
	TIndividual();
	virtual ~TIndividual();

	vector <int> x_var;
	vector <double> y_obj;

	void   rnd_init();
	void   obj_eval();

    bool   operator<(const TIndividual &ind2);
    bool   operator==(const TIndividual &ind2);
    void   operator=(const TIndividual &ind2);

	void show_objective();
	void show_variable();

};


#endif


