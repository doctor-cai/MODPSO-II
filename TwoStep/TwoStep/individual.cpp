#include "individual.h"

TIndividual::TIndividual()
{
	for(int i=0; i<numVariables; i++)
		x_var.push_back(i+1);	//this is better
	//	x_var.push_back(rnd_uni(&rnd_uni_init) * numVariables);

	for(int n=0; n<numObjectives; n++)
        y_obj.push_back(0.0);
}

TIndividual::~TIndividual()
{

}

void TIndividual::rnd_init()
{	
	//*初始化粒子位置
	IGLP(x_var,node);

//	heuristic_initial(x_var,node);
//	random_init(x_var);
//	rehearse(x_var);

}

void TIndividual::obj_eval()
{
    objectives(x_var,y_obj);
}

void TIndividual::operator=(const TIndividual &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
}

void TIndividual::show_objective()
{
    for(int n=0; n<numObjectives; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void TIndividual::show_variable()
{
    for(int n=0; n<numVariables; n++)
		printf("%d ",x_var[n]);
	printf("\n");
}

bool TIndividual::operator<(const TIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<numObjectives; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false; //min
	//	if(ind2.y_obj[n]>y_obj[n]) return true;  //max
	}
	if(ind2.y_obj==y_obj) return false;
//	if(ind2.y_obj==y_obj) return true;
	return dominated;
}


bool TIndividual::operator==(const TIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}
