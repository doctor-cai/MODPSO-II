
#include "moead.h"

TMOEAD::TMOEAD()
{
	velocity = new double *[popsize];
	for (int i = 0; i < popsize; i++)
	{
		velocity[i] = new double[numVariables];
	}

	idealpoint = new double[numObjectives];
	indivpoint = new TIndividual[numObjectives];    
	// initialize ideal point	
    for (int n=0; n<numObjectives; n++) 
	{
		if (optimization == 0)
		{
			idealpoint[n] = 1.0e+30;   // for min
		} 
		else if (optimization == 1)
		{
			idealpoint[n] = -1.0e+30;  // for max
		}
	
		indivpoint[n].rnd_init();
		indivpoint[n].obj_eval();
	}
}

TMOEAD::~TMOEAD()
{
	for (int i=0; i<popsize; i++)
	{
		delete []velocity[i];
	}
	delete []velocity;

	delete []idealpoint;    
	delete []indivpoint;
}

// indiv stands for the X(the number of the upindex is N),
//rnd_init() is used to generate the x(the number of the downindex equals to the number of the objectives)

void TMOEAD::init_population()
{
    for (int i=0; i<popsize; i++)
	{
		population[i].indiv.rnd_init();
		population[i].indiv.obj_eval();
		update_reference(population[i].indiv);
	}
}

void TMOEAD::init_velocity()
{
    for (int i=0; i<popsize; i++)
	{
		for (int j=0; j<numVariables; j++)
		{
			velocity[i][j] = 0.0;		
		}
	}
}

// Store personal bests (both variable and fitness values) of particles

void TMOEAD::store_pbests()  
{
	unsigned int i;
	
	// Store variable values of personal bests 
	for (i=0; i<popsize; i++)
	{
		pbest.push_back(population[i].indiv);
	}		
}

// Compute new velocity of each particle in the population 

void TMOEAD::compute_velocity(int index, TIndividual &child, int gen) 
{
	unsigned int i,j, k, m, l;

	child = population[index].indiv;
	TIndividual NeighborBest;

	/***************版本1 从邻居随机选择作为gbest**********************/

    int RandNum = rnd_uni(&rnd_uni_init) * population[index].table.size();
	int RandNeighbor = population[index].table[RandNum];
	NeighborBest = population[RandNeighbor].indiv;
	
	/*********************version 2 从邻居选择加权最小的作为gbest***********/

/*	double minaggre = 1000000;
	int    minaggre_index = 0;
	for (i = 0; i < population[index].table.size(); i++)
	{
		int index_n = population[index].table[i];
		double aggregation = population[index_n].indiv.y_obj[0] * population[index_n].namda[0] +
							 population[index_n].indiv.y_obj[1] * population[index_n].namda[1];
		if (minaggre > aggregation)
		{
			minaggre = aggregation;
			minaggre_index = i;
		}
	}

	NeighborBest = population[minaggre_index].indiv;*/

	/********************* version 3 从邻居中选择非支配的邻居作为gbest****************/
/*	int index_n[10];
	int *num_dominated = new int[popsize];

	for (i = 0; i < population[index].table.size(); i++)
	{
		index_n[i] = population[index].table[i];
	}
	
	
    for (i=0; i<popsize; i++)
	{
		num_dominated[i] = 0;
	}

	for (i = 0; i < population[index].table.size(); i++)
	{
		for (j = i+1; j < population[index].table.size(); j++)
		{
			int temp = donimate_judge(population[index_n[i]].indiv.y_obj, population[index_n[j]].indiv.y_obj);
			if (temp == -1) //i donimate j
				 num_dominated[j]++;
			 else if (temp == 1)//j donimate i
				 num_dominated[i]++;
		}

		int random1=index_n[i];
		if (num_dominated[i] == 0)
		{
			NeighborBest = population[random1].indiv;
		}
	}*/
	

	/************* version 4 sigma method to chose the gbest from the neighborhoods***/

//  shit! 为什么不能用vector 容器存储，一用就错！！！
//	vector<double> deta_obj;
//	vector<double> deta_distance;
	
/*	double deta_obj[100];
	double deta_distance[100];

	for (i=0; i<population[index].table.size(); i++)
	{
		int tab = population[index].table[i];
		double yzhi = population[tab].indiv.y_obj[0];
		
		double a = population[tab].indiv.y_obj[0]*population[tab].indiv.y_obj[0];
		double b = population[tab].indiv.y_obj[1]*population[tab].indiv.y_obj[1];
		double c = a - b;
		double d = a + b;
		double e = c/d;
		deta_obj[i] = c/d;
	}

	double a1 = population[index].indiv.y_obj[0]*population[index].indiv.y_obj[0];
	double b1 = population[index].indiv.y_obj[1]*population[index].indiv.y_obj[1];
	double deta_obj_index = (a1 - b1)/(a1 + b1);	

	for (i=0; i<population[index].table.size(); i++)
	{
		double yt = fabs(deta_obj[i] - deta_obj_index);
		deta_distance[i] = yt;
	}

	double temp = deta_distance[0];
	int    distance_index = 0;

	for (i=1; i<population[index].table.size(); i++)
	{
		if (deta_distance[i] <= temp)
		{
			temp = deta_distance[i];
			distance_index = i;
		}
	}

	int e = population[index].table[distance_index];
	
	for (i=0; i<numVariables; i++)
	{
		NeighborBest.x_var[i] = population[e].indiv.x_var[i];
	}

	for (i=0; i<numObjectives; i++)
	{
		NeighborBest.y_obj[i] = population[e].indiv.y_obj[i];
	}*/

	//NeighborBest = population[e].indiv;

	for(j=0; j<numVariables; j++)
	{
		int v1, v2;

		if (pbest[index].x_var[j] == population[index].indiv.x_var[j]) 
		{
			v1 = 0;
		}
		else 
		{
			v1 = 1;
		}

		if (NeighborBest.x_var[j] == population[index].indiv.x_var[j]) 
		{
			v2 = 0;
		}
		else
		{
			v2 = 1;
		}

// 		velocity[i][j] = 0.4 * velocity[i][j] + 1.0 * rnd_uni(&rnd_uni_init) * 
//  		v1 + 1.0 * rnd_uni(&rnd_uni_init) * v2;

// 		velocity[i][j] = 0.4 * velocity[i][j] + 0.5 * rnd_uni(&rnd_uni_init) * 
//  		v1 + 0.5 * rnd_uni(&rnd_uni_init) * v2;

		velocity[index][j] = rnd_uni(&rnd_uni_init) * velocity[index][j] + 
							1.494 * rnd_uni(&rnd_uni_init) * v1 + 
							1.494 * rnd_uni(&rnd_uni_init) * v2;

// 		velocity[index][j] = (0.329 * (99 - gen)/99 + 0.4) * velocity[index][j] + 
// 			1.494 * rnd_uni(&rnd_uni_init) * v1 + 
// 							1.494 * rnd_uni(&rnd_uni_init) * v2;
//      this is inferior to the above one
		

		double sigmoid = 1.0/(1.0+exp(-velocity[index][j]));
		
		if (rnd_uni(&rnd_uni_init) < sigmoid)
		{
			velocity[index][j] = 1;
		}
 		else 
		{
			velocity[index][j] = 0;
		}
	}

  /* Calculate new positions of particles */
//       popVar[i][j] = popVar[i][j] + velocity[i][j];	

///////////////////////////////////////////////////////////////////////////////////////////////
// select the label indentifier of a node which has largest degree
//   for(i = 0; i < popsize; i++){
//     for(j = 0; j < numVariables; j++){
// 		if (velocity[i][j] == 1){
// 			int neighborsize = node[j].neighbours.size();
// 			
// 			if(neighborsize > 1){
// 				vector<int> neighbor_degree;		//存储i节点邻居
// 				for (int nei = 0; nei < neighborsize; nei++){
// 					int neighbornode = node[j].neighbours[nei];
// 					neighbor_degree.push_back(node[neighbornode].degree);
// 				}
// 
// 				int max_degree = *max_element(neighbor_degree.begin(),neighbor_degree.end());
// 				int nodeindex;
// 				
// 				for (int tempindex = 0; tempindex < neighborsize; tempindex++){
// 					if (max_degree == neighbor_degree[tempindex]){
// 						nodeindex = tempindex;
// 						break;
// 					}
// 				}
// 				
// 				popVar[i][j] = popVar[i][node[j].neighbours[nodeindex]];
// 				
// 			}//end if neighbor.size > 1
// 			else {
// 				if (neighborsize == 1){
// 					popVar[i][node[j].neighbours[0]] = popVar[i][j];
// 				} 
// 				else{
// 					//	cout<<"第"<<i+1<<"个节点无连接"<<endl;
// 					//	pos[i] = 0;
// 				}
// 				
// 			}
// 		}//end if
// 	}//end i
//   }//end j
///////////////////////////////////////////////////////////////////////////////////////////////

/*********************************************************************************************/
// select the dominated label indentifier 

	  for (j=0; j<numVariables; j++)
	  {
		  if (velocity[index][j] == 1)
		  {
				int	sum = 0;
				int maxr = -1;//record index of i's neighbour which ...	
				int label = -1;
				int	temp = 1;
							
				if (node[j].neighbours.size() > 1)
				{
					for (m=0; m<node[j].neighbours.size(); m++)
					{
						int counter = 1;//record no. of nodes that has same label with j
									
						for (k=m+1; k<node[j].neighbours.size(); k++ )
						{
							int p = child.x_var[node[j].neighbours[m]];
							int q = child.x_var[node[j].neighbours[k]];
							if( p == q )  
							{
								counter++;
							} 
						}//end k
									
						if (temp < counter)
						{
							maxr = m;
							temp = counter;
						}
					}//end m
								
					for (l=0; l<node[j].neighbours.size(); l++)
					{
						int u = child.x_var[node[j].neighbours[l]];
						int v = child.x_var[j];

						if ( u == v )	
						{
							label = u;
						}
					}//endl

					if (label != -1 && maxr == -1)
					{
						child.x_var[j] = label;
					}
					else 
					{
						if (maxr != -1)
						{
							child.x_var[j] = child.x_var[node[j].neighbours[maxr]];
						} 
						else
						{
							double r3 = rnd_uni(&rnd_uni_init);
							int randneighbor = r3*node[j].neighbours.size();

							if (r3 == 1)
							{
								child.x_var[j] = child.x_var[node[j].neighbours[randneighbor-1]];
							}
							else
							{
								child.x_var[j] = child.x_var[node[j].neighbours[randneighbor]];				
							}
						}			
					}

				}//end if node[i].neighbours.size() > 1
				else 
				{
					if (node[j].neighbours.size() == 1)
					{
						child.x_var[j] =  child.x_var[node[j].neighbours[0]];	// this is better
					//	child.x_var[node[j].neighbours[0]] =  child.x_var[j];
					} 
					else
					{
			//			cout<<"the "<<i+1<<"th node has no connections"<<endl;
	//					pos[i] = 0;
					}
				}				
		  }//end if velocity[i][j] == 1	 
	}//end j
}//end function

void TMOEAD::update_pbests(int index) /* Update personal bests of particles in the population */
{
	unsigned int j, sum, better, counter;
	
	sum = 0; 
	counter = 0;
	better = 0;

	for(j = 0; j < numObjectives; j++)
	{
		if( ((population[index].indiv.y_obj[j] <= pbest[index].y_obj[j]) && (optimization == 0)) 
			|| ((population[index].indiv.y_obj[j] >= pbest[index].y_obj[j]) && (optimization == 1)))
		{
			sum += 1;
		}

		if( ((population[index].indiv.y_obj[j] < pbest[index].y_obj[j]) && (optimization == 0)) 
			|| ((population[index].indiv.y_obj[j] > pbest[index].y_obj[j]) && (optimization == 1)))
		{
			counter += 1;
		}
	}	
		
	if (sum == numObjectives) 		
	{ /* current pop dominates pbest */
		better = 0;
	} 
	else 
	{
		if (sum == 0)           /* pbest dominates current pop */
		{
			better = 1;
		}
		else if (counter == 1)
		{
			double temp1 = population[index].namda[0]*population[index].indiv.y_obj[0]+
							population[index].namda[1]*population[index].indiv.y_obj[1];
			double temp2 = population[index].namda[0]*pbest[index].y_obj[0]+
							population[index].namda[1]*pbest[index].y_obj[1];
			
			if (((temp1<temp2) && (optimization == 0)) || ((temp1>temp2) && (optimization == 1)))
			{
				better = 0;//*/
			}
		}
		//	better = rnd_uni(&rnd_uni_init); /* both are nondominated,randomly select one */
	}
		
	if (better == 0)
	{
		//pbest[index] = population[index].indiv;
		for(j = 0; j < numObjectives; j++)
		{
			pbest[index].y_obj[j] = population[index].indiv.y_obj[j];
		}

		for(j = 0; j < numVariables; j++)
		{
			pbest[index].x_var[j] = population[index].indiv.x_var[j];				
		}
	}
}


void TMOEAD::update_pbests_child(TIndividual &ind, int index)
{
	unsigned int j, sum, better, counter;
	
	sum = 0; counter = 0;

	for(j=0; j<numObjectives; j++)
	{
		if( ((ind.y_obj[j] <= pbest[index].y_obj[j]) && (optimization == 0)) 
			|| ((ind.y_obj[j] >= pbest[index].y_obj[j]) && (optimization == 1)))
		{
			sum += 1;
		}

		if( ((ind.y_obj[j] < pbest[index].y_obj[j]) && (optimization == 0)) 
			|| ((ind.y_obj[j] > pbest[index].y_obj[j]) && (optimization == 1)))
		{
			counter += 1;
		}
	}	
		
	if (sum == numObjectives) 		
	{ /* current pop dominates pbest */
		better = 0;
	} 
	else 
	{
		if (sum == 0)           /* pbest dominates current pop */
		{
			better = 1;
		}
		else if (counter == 1)
		{
			//*
			if (rnd_uni(&rnd_uni_init) < 0.5)
			{
				better = 1;
			}
			else 
			{
				better = 0;//*/
			}

			/*
			double temp1 = population[index].namda[0]*ind.y_obj[0]+
							population[index].namda[1]*ind.y_obj[1];
			double temp2 = population[index].namda[0]*pbest[index].y_obj[0]+
							population[index].namda[1]*pbest[index].y_obj[1];
			if (((temp1<temp2) && (optimization == 0)) || ((temp1>temp2) && (optimization == 1)))
				better = 0;//*/
		}
	}
		
		if (better == 0)
		{
//			pbest[index] = population[index].indiv;
			for(j = 0; j < numObjectives; j++)
				pbest[index].y_obj[j] = ind.y_obj[j];
			for(j = 0; j < numVariables; j++)
 				pbest[index].x_var[j] = ind.x_var[j];				
		} 
}

// initialize a set of evenly-distributed weight vectors 
/***************version 1 MODPSO.weightvectors*******************/

void TMOEAD::init_uniformweight()
{   
    for(int i=0; i<popsize; i++)
	{
        TSOP sop;		    //when definite the sop class, erase the sop's initial value
		sop.array.push_back(i);
		sop.array.push_back(popsize-i);
		for(int j=0; j<sop.array.size(); j++)
			sop.namda.push_back(1.0*sop.array[j]/popsize);  
		population.push_back(sop); 		             
	}
}

/***************version 2 MODPSO.weightvectors*******************/

/*void TMOEAD::init_uniformweight(int sd)
{   
    for(int i=0; i<=sd; i++)
	{
        TSOP sop;
		double a;
		a=1.0*i/sd;
		sop.namda.push_back(a);
		sop.namda.push_back(1-a);
		population.push_back(sop); 		
	}
	pops = population.size();
}*/

void TMOEAD::update_weightvector()
{
	double sum=0.0;
	for (int i=0; i<popsize; i++)
	{
		double acc=population[i].namda[0];
		for (int j=1; j<numObjectives; j++)
		{
			acc=acc * population[i].namda[j];
		}
		if (acc!=0)
		{
			for (int j=0; j<numObjectives; j++)
			{
				
				sum = sum + population[i].namda[j];
			}
			for (int j=0; j<numObjectives; j++)
			{
				population[i].namda[j] = 1.0/population[i].namda[j];
				population[i].namda[j] = population[i].namda[j]/sum;
			}
		}
	}
}

// initialize the neighborhood of subproblems based on the distances of weight vectors
// choose the inche neighborhoods based on the Euclidean distances of weight vectors 

void TMOEAD::init_neighbourhood()
{
    double *x   = new double[popsize];
	int    *idx = new int[popsize];
	for(int i=0; i<popsize; i++)
	{	
		for(int j=0; j<popsize; j++)
		{
		    x[j]    = distanceVector(population[i].namda,population[j].namda);
			idx[j]  = j;			
		}
		minfastsort(x,idx,popsize,niche);   
		for(int k=0; k<niche; k++)   
		{
			//table: store the neighbourhoods' index
			population[i].table.push_back(idx[k]);  
		}
	/*	minfastsort(x,idx,pops,Tr);   
		for(  k=0; k<Tr; k++)
		{
			//table: store the neighbourhoods' index
			population[i].tabler.push_back(idx[k]);  
		}

*/
	}
    delete [] x;
	delete [] idx;

}



//void TMOEAD::update_problem(TIndividual &indiv, int id,int Tr)
/*void TMOEAD::update_problem(TIndividual &indiv, int id)
{
	//新解已产生，新的函数值也计算出，求新解Xnew对应的各个粒子的函数值，
	//取最小函数值对应的子问题的index;
	// id 怎么用
	double *f1 = new double[popsize];
	//double temp;
	int m=0;
	int k;

	for ( k=0; k<popsize; k++)
	{
		f1[k]=scalar_func(indiv.y_obj, population[k].namda, indivpoint);
		
	}

	//      step1  求最小的值对应的index 
	
	//temp=f1[0];
	for ( k=1; k<popsize; k++)
	{
		if (f1[m]>f1[k])
		{
			m=k;      //m中存放距离新解Xnew近的子问题的index
		}
	}
	//  step2 确定replacement neighborhood Tr 子问题i的Tr个closest subproblems
	//  初始化邻居时，已经把每个子问题的最近邻居求出 只不过个数为niche 现在要求是Tr个
	double *x   = new double[popsize];
	int    *idx = new int[popsize];
	for(int i=0; i<popsize; i++)
	{	
		for(int j=0; j<popsize; j++)
		{
		    x[j]    = distanceVector(population[i].namda,population[j].namda);
			idx[j]  = j;			
		}
		minfastsort(x,idx,popsize,Tr);   
		for( k=0; k<Tr; k++)   
			population[i].tabler.push_back(idx[k]);  //table: store the neighbourhoods' positions存储邻居的index

	}

  // step3  Tr个邻居被更新
	for ( k=0;k<Tr;k++)
	{
		int    m1  = population[m].tabler[k];
		double f11, f12;
		f11 = scalar_func(population[m1].indiv.y_obj, population[m1].namda, indivpoint);
		f12 = scalar_func(indiv.y_obj, population[m1].namda, indivpoint);

		if(f12<f11) population[m1].indiv = indiv;
	}
	delete [] f1;
	delete [] x;
	delete [] idx;
}*/


/*void TMOEAD::update_problem(TIndividual &indiv, int id)
{
    for(int i=0; i<niche; i++)
	{
		int    k  = population[id].table[i];
		double f1, f2;
		f1 = scalar_func(population[k].indiv.y_obj, population[k].namda, indivpoint);
		f2 = scalar_func(indiv.y_obj, population[k].namda, indivpoint);

		if(f2<f1) population[k].indiv = indiv;		
	}
}*/

//  MOEA/D DE
void TMOEAD::update_problem(TIndividual &indiv, int id)
{
	int    j,k;
	int    neigh1;
	int    neigh2;
	int    c      = 0;

//	vector <double> x1;
//	vector <int> idx1;
//	vector <double> x2;
//	vector <int> idx2;
	
	double *f11   = new double[niche];
	double *f12   = new double[niche];
	double *f21   = new double[popsize];
	double *f22   = new double[popsize];

	double *x1    = new double[niche];
	int    *idx1  = new int[niche];
	double *x2    = new double[popsize];
	int    *idx2  = new int[popsize];

	double random = rnd_uni(&rnd_uni_init);

	if (random < deta)  //p取邻居
	{
		for (j=0; j<niche; j++)
		{
			neigh1  = population[id].table[j];
			f11[j]  = scalar_func(indiv.y_obj, population[neigh1].namda, indivpoint);
			f12[j]  = scalar_func(population[neigh1].indiv.y_obj, population[neigh1].namda, indivpoint);
			if (f11[j] <= f12[j])
			{
				c = c + 1;
			}
		}
		if (c <= nr)
		{
			for (j=0; j<niche; j++)
			{
				if (f11[j] <= f12[j])
				{
					population[population[id].table[j]].indiv = indiv;
				}
			}
		}
		else
		{
			int index1 = 0;

			for (j=0; j<niche; j++)
			{
				
				if (f11[j] <= f12[j])
				{
					neigh2  = population[id].table[j];  
					x1[index1] = distanceVector(population[neigh2].indiv.y_obj,indiv.y_obj);
					idx1[index1] = j;
					index1++;
				}
			}

			minfastsort(x1, idx1, index1, nr);
			for (k=0; k<nr; k++)   
			{
				population[id].tabler.push_back(idx1[k]);
				population[population[id].tabler[k]].indiv = indiv;
			}
			
		}	
	}//end if
	else
	{
		c = 0;
		for (j=0; j<popsize; j++)
		{
			f21[j] = scalar_func(indiv.y_obj, population[j].namda, indivpoint);
			f22[j] = scalar_func(population[j].indiv.y_obj, population[j].namda, indivpoint);
			if (f21[j] <= f22[j])
			{
				c = c + 1;
			}
		}
		if (c <= nr)
		{
			for (j=0; j<popsize; j++)
			{
				if (f21[j] <= f22[j])
				{
					population[j].indiv = indiv;
				}
			}
		}
		else
		{
			int index2=0;
			for (j=0; j<popsize; j++)
			{
				if (f21[j] <= f22[j])
				{
					x2[index2] = distanceVector(population[j].indiv.y_obj, indiv.y_obj);
					idx2[index2] = j;
					index2++;
				}
			}
			minfastsort(x2, idx2, index2, nr);
			for (k=0; k<nr; k++)   
			{
//				population[id].tabler.push_back(idx2[k]);
//				population[population[id].tabler[k]].indiv = indiv;
				population[idx2[k]].indiv = indiv;
			}
		}
	}
//	delete [] f1;
	delete []f11;
	delete []f12;    
	delete []f21;
	delete []f22;
//	delete []x1;
//	delete []idx1;
//	delete []x2;
//	delete []idx2;
}




/*void TMOEAD::update_problem(TIndividual &indiv, int id)
{
	
	//compute the distance between the new objectives value and other particles

	double *x1   = new double[popsize];
	int k;
//	int    *idx1 = new int[pops];
	for (int i=0;i<popsize;i++)
	{
			x1[i]=distanceVector(population[i].indiv.y_obj,indiv.y_obj);
		//	idx1[j]=j;
	}

	int i1=0;
	for (i=1;i<popsize;i++)
	{
		if (x1[i1]>x1[i])
		{
			i1=i;
		}
	}

	//  step2 compute replacement neighborhood Tr, subproblem i's Tr closest subproblems
	double *x   = new double[popsize];
	int    *idx = new int[popsize];
	for(int q=0; q<popsize; q++)
	{	
		for(int j=0; j<popsize; j++)
		{
		    x[j]    = distanceVector(population[q].namda,population[j].namda);
			idx[j]  = j;			
		}

		minfastsort(x,idx,popsize,Tr);
		
		for( k=0; k<Tr; k++)  
		{
			//table: store the neighbourhoods' index
			population[q].tabler.push_back(idx[k]);  
		}

	}
    

  // step3  update the Tr subproblems
	for (k=0;k<Tr;k++)
	{
		//把i写成id,结果还可以
		int    m  = population[i1].tabler[k];
		double f11, f12;
		f11 = scalar_func(population[m].indiv.y_obj, population[m].namda, indivpoint);
		f12 = scalar_func(indiv.y_obj, population[m].namda, indivpoint);

		if(f12<f11) population[m].indiv = indiv;
	}
	delete [] x1;
	delete [] x;
	delete [] idx;

//	delete [] idx1;
}*/

/*void TMOEAD::update_problem(TIndividual &indiv, int id)
{
	TIndividual indiv1;
	double m1=rnd_uni(&rnd_uni_init)*population[id].table.size();
	double m2=rnd_uni(&rnd_uni_init)*population[id].table.size();
	int m3,m4;
	m1=int (m1);
	m2=int (m2);
	int k;
	double *f1=new double[pops];
	if (m1!=m2)
	{
		m3=population[id].table[m1];
		m4=population[id].table[m2];

		//generate the Xnew
	//	realbinarycrossover1 (population[m3].indiv, population[m4].indiv, indiv);
		indiv.obj_eval();  //indiv.y_obj[i]已经计算出
	}

	//计算新解的目标函数值和其他粒子的目标函数值的距离

	double *x1   = new double[pops];
//	int    *idx1 = new int[pops];
	for (int i=0;i<pops;i++)
	{
			x1[i]=distanceVector(population[i].indiv.y_obj,indiv.y_obj);
		//	idx1[j]=j;
	}
	double x11=x1[0];
	int i1=0;
	for (i=1;i<pops;i++)
	{
		if (x11>x1[i])
		{
			x11=x1[i];
			i1=i;
		}
	}
	for( k=0; k<pops; k++)
	{
		f1[k]=scalar_func(population[k].indiv.y_obj, population[k].namda, indivpoint);
	}
	
	//      step1  求最小的值对应的index 
	double temp;
	int a1=0;
	temp=f1[0];
	for (k=1; k<pops; k++)
	{
		if (temp>f1[k])
		{
			temp=f1[k];
			a1=k;
		}
	}
	//  step2 compute replacement neighborhood Tr, subproblem i's Tr closest subproblems
	double *x   = new double[pops];
	int    *idx = new int[pops];
	for(int q=0; q<pops; q++)
	{	
		for(int j=0; j<pops; j++)
		{
		    x[j]    = distanceVector(population[q].namda,population[j].namda);
			idx[j]  = j;			
		}
		minfastsort(x,idx,pops,Tr);   
		for( k=0; k<Tr; k++)   
			population[q].tabler.push_back(idx[k]);  //table: store the neighbourhoods' index

	}
    

  // step3  update the Tr subproblems
	for (k=0;k<Tr;k++)
	{
		//把i写成id,结果还可以
		int    m5  = population[i1].tabler[k];
		double f11, f12;
		f11 = scalar_func(population[m5].indiv.y_obj, population[m5].namda, indivpoint);
		f12 = scalar_func(indiv.y_obj, population[m5].namda, indivpoint);

		if(f12<f11) population[m5].indiv = indiv;
	}
	delete [] f1;
	delete [] x;
	delete [] idx;
	delete [] x1;
//	delete [] idx1;
}*/
// update the reference point
void TMOEAD::update_reference(TIndividual &ind)
{
	for(int n=0; n<numObjectives; n++)    
	{
		if (optimization == 0)
		{
			if(ind.y_obj[n] < idealpoint[n])
			{
				idealpoint[n]  = ind.y_obj[n];
				indivpoint[n]  = ind;
			}
		} 
		else if (optimization == 1)
		{
			if(ind.y_obj[n] > idealpoint[n])
			{
				idealpoint[n]  = ind.y_obj[n];
				indivpoint[n]  = ind;
			} 
		}		
	}
}

void TMOEAD::evolution(int gen)
{
	
    for(int i=0; i<popsize; i++)
	{			
		int   n  =  i; 		
		TIndividual child, child2;

		compute_velocity(n,child,gen);

		if (i < popsize*mutate_posibility)  
//		realmutation1(child); //师兄原版程序扰动
		realmutation(child);
		child.obj_eval();
		update_problem(child, n);
		update_reference(child);
		update_pbests(n);
	}

}


void TMOEAD::run(int mg, int rn, int ithData, int ithRun)
{
    // sd: integer number for generating weight vectors
	// nc: size of neighborhood
	// mg: maximal number of generations 
	
	init_uniformweight();
	
	//update_weightvector();
    init_neighbourhood();

	init_population();
	
	store_pbests();
	init_velocity();

	for(int gen = 2; gen <= mg; gen++)   
	{
	//	cout<<"the "<<gen<<"th iteration finished"<<endl;
		evolution(gen);
	}

	char savefilename0[1024];
	char savefilename1[1024];

	if (SignedFlag == 0)
	{
		sprintf(savefilename0, "PF_unsigned/TwoStep_data%d_PF%d.txt", ithData, ithRun);
		sprintf(savefilename1, "variable_unsigned/TwoStep_data%d_X%d.txt", ithData, ithRun);
	}

	if (SignedFlag == 1)
	{
		sprintf(savefilename0, "PF_signed/TwoStep_data%d_PF%d.txt", ithData, ithRun);
		sprintf(savefilename1, "variable_signed/TwoStep_data%d_X%d.txt", ithData, ithRun);
	}

	save_results(savefilename0, savefilename1);

	population.clear();
}

void TMOEAD::save_results(char savePF[100], char saveChrom[100])
{	
	vector<vector<double> >pop_y;
	vector<int> index;
	vector<vector<double> >front;
	int i,j,k;

	ofstream fout_front(savePF, std::ios::out);
	ofstream fout_chrom(saveChrom, std::ios::out);

	for (i=0; i<popsize; i++)
	{
		pop_y.push_back(population[i].indiv.y_obj);
	}

	GetNondominatedSolution(pop_y, front, index) ;
	
	for (int n = 0; n<front.size(); n++)
	{
		for (int k = 0; k<numObjectives; k++)
			fout_front << front[n][k] << "  ";
		fout_front << endl;
	}
	fout_front.close();


	for (k=0; k<index.size(); k++)
	{
		int Anum=0;
		int countn = 0, countp = 0;

		cout << "pop individual " << setw(3) << left << index[k] << " -- ";
		fout_chrom << "pop individual " << setw(3) << left << index[k] << " -- ";
		double Q = calcQ(population[index[k]].indiv.x_var, Anum, countn, countp);
		
		if (check_label())
		{
			double temp_nmi = calc_NMI(population[index[k]].indiv.x_var, LabelName);
			fout_chrom << " NMI = " << temp_nmi << setprecision(5) << "    ";
			cout << " NMI = " << temp_nmi << setprecision(5) << "    ";
			NMI.push_back(temp_nmi);
		}

		cout << "Q = " << Q << setprecision(5) << " ";
		cout << "countn = " << setw(7) << left << countn << " ";
		cout << "countp = " << setw(7) << left << countp << " ";
		fout_chrom << "Q = " << Q << setprecision(5) << " ";
		fout_chrom << "countn = " << setw(7) << left << countn << " ";
		fout_chrom << "countp = " << setw(7) << left << countp << " ";
		cout << "cluster is " << setw(7) << left << Anum << " ";
		fout_chrom << "cluster is " << setw(7) << left << Anum << " ";
		

		modularity.push_back(Q);

		int energyfunc=0,m;

//		computing the energyfunction
		for (i=0; i<numVariables; i++)
		{
			for (j=i+1; j<numVariables; j++)
			{
				if (population[index[k]].indiv.x_var[i] == population[index[k]].indiv.x_var[j])
				{
					m = 1;
				}
				else 
				{
					m = -1;
				}
//				energyfunc += (1 - AdjacentMatrix[i][j]*m)/2;

				int aij = 0;
				if (AdjacentMatrix[i][j]==0x31)
				{
					aij = 1;
				}
				else 
				{
					if (AdjacentMatrix[i][j]==0x2d)	
						aij = -1;
				}
				energyfunc += (1 - aij*m)/2;
			}
		}

		cout << "hs = " << setw(7) << left << energyfunc << endl;
		fout_chrom << "hs = " << setw(7) << left << energyfunc << endl;

		fout_chrom << "variables is" << endl;
		for (i = 0; i<numVariables; i++)
		{
			fout_chrom << population[index[k]].indiv.x_var[i] << " ";
		}
		fout_chrom << endl << endl;

	}
	if (check_label())
	{
	//	cout << "max NMI = " << *max_element(NMI.begin(), NMI.end()) << "    ";
	//	fout_chrom << "max NMI = " << *max_element(NMI.begin(), NMI.end()) << endl;
		pm_NMI.push_back(*max_element(NMI.begin(), NMI.end()));
		Thirty_Run_NMI.push_back(*max_element(NMI.begin(), NMI.end()));
	}

//	cout << "max modularity = " << *max_element(modularity.begin(), modularity.end()) << endl;
//	fout_chrom << "max modularity = " << *max_element(modularity.begin(), modularity.end()) << endl << endl;
	pm_modularity.push_back(*max_element(modularity.begin(), modularity.end()));
	Thirty_Run_modularity.push_back(*max_element(modularity.begin(), modularity.end()));

	fout_front.close();

	pop_y.clear();
	index.clear();
	front.clear();

}



void TMOEAD::operator=(const TMOEAD &emo)
{
	population  = emo.population;
	indivpoint  = emo.indivpoint;
//	niche       = emo.niche;
} 

