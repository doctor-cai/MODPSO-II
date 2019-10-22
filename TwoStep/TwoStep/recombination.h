#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "individual.h"


template <class T>  //execute two way crossover
void realbinarycrossover(T &parent1, T &parent2, T& child1, T& child2)
{
	for (short i = 0; i < numVariables; i++)
	{
		child1.x_var[i] = parent1.x_var[i];
		child2.x_var[i] = parent2.x_var[i];
	}

	if (rnd_uni(&rnd_uni_init) <= cross_posibility)
	{

		short position1 = short(rnd_uni(&rnd_uni_init)*numVariables);
		short position2 = short(rnd_uni(&rnd_uni_init)*numVariables);
		short a = child1.x_var[position1];
		short b = child2.x_var[position2];

		vector<short>indexA; //store index in cha that has same gene with position
		vector<short>indexB; //store index in chb that has same gene with position

		for (short m = 0; m < numVariables; m++)
		{
			if (child1.x_var[m] == a)
			{
				indexA.push_back(m);
			}
			if (child2.x_var[m] == b)
			{
				indexB.push_back(m);
			}
		}

		for (short j = 0; j < indexA.size(); j++)
		{
			child2.x_var[indexA[j]] = a;
		}

		for (short n = 0; n < indexB.size(); n++)
		{
			child1.x_var[indexB[n]] = b;
		}
	}
}

// generate a random between the [-1,1],update the individual position
template <class T>
void rand_update_indivi(T &child)
{
	for (short i = 0; i<numVariables; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= 0.01)
		{
			double random1 = rnd_uni(&rnd_uni_init);
			double random2 = -rnd_uni(&rnd_uni_init);
			if (rnd_uni(&rnd_uni_init) <= 0.5)
			{
				child.x_var[i] = child.x_var[i] + random1*child.x_var[i];
			}
			else
			{
				child.x_var[i] = child.x_var[i] + random2*child.x_var[i];
			}
			//double random1=rnd_uni(&rnd_uni_init);
			//child.x_var[i]=child.x_var[i]+random1*child.x_var[i];
		}
	}
}

template <class T>  //execute two way crossover
void realbinarycrossover1(T &parent1, T &parent2, T& child)
{
	T child1, child2;
	for (short i = 0; i < numVariables; i++)
	{
		child1.x_var[i] = parent1.x_var[i];
		child2.x_var[i] = parent2.x_var[i];
	}
	if (rnd_uni(&rnd_uni_init) <= cross_posibility)
	{
		for (short j = 0; j < numVariables; j++)
		{
			double random1 = rnd_uni(&rnd_uni_init);
			double random2;
			random2 = random1*child1.x_var[j] + (1 - random1)*child2.x_var[j];
			child.x_var[j] = short(random2);
		}
	}
	/*if (rnd_uni(&rnd_uni_init) <= cross_posibility)
	{

	short position1 = short(rnd_uni(&rnd_uni_init)*numVariables);
	short position2 = short(rnd_uni(&rnd_uni_init)*numVariables);
	short a = child1.x_var[position1];
	short b = child2.x_var[position2];

	vector<short>indexA; //store index in cha that has same gene with position
	vector<short>indexB; //store index in chb that has same gene with position

	for (short m = 0;m < numVariables;m++)
	{
	if (child1.x_var[m] == a )
	{
	indexA.push_back(m);
	}
	if ( child2.x_var[m] == b )
	{
	indexB.push_back(m);
	}
	}

	for ( short j = 0;j < indexA.size();j++ )
	{
	child2.x_var[indexA[j]] = a;
	}

	for ( short n = 0;n < indexB.size();n++ )
	{
	child1.x_var[indexB[n]] = b;
	}
	/*double aggregation1 = child1.y_obj[0] *child1.namda[0] +
	child1.indiv.y_obj[1] * child1.namda[1];
	double aggregation2 = child2.y_obj[0] *child2.namda[0] +
	child2.indiv.y_obj[1] * child2.namda[1];

	}*/
}


template <class T>
void realmutation(T &ind)  //NBM mutation
{
	for (short j = 0; j < numVariables; j++)
	{
		if (rnd_uni(&rnd_uni_init) <= mutate_posibility)
		{
			short indentifier = ind.x_var[j];

			short rand_v_ns = node[j].neighbours.size();
			short negative_n = node[j].neighbours_n.size();

			//版本1,只针对positive edges,node 的邻居,
			for (short i = 0; i < rand_v_ns; i++)
			{
				short neighborX = node[j].neighbours[i];
				ind.x_var[neighborX] = indentifier;
			}

		}
	}

	return;
}



/********************************************************/
/* local search procedure:hill climbing process(greedy) */
/* strictly follow the version in the paper:meme-net.pdf*/
/* randomly select a chrom and a gene,say the chosen    */
/* chrom has m clusters,and the chosen gene belongs to  */
/* the n-th cluster,assign the gene to other clusters,  */
/* thus we get m-1 new clusters.if the best D value in  */
/* the m-1 clusters is larger than the chosen chrom's   */
/* copy the new chrom to the chosen chrom				*/
/********************************************************/
//*
template <class T>
void LocalSeach(T &ind, vector<double>&namda)
{
	short i, j, k;
	short RandGene = numVariables * rnd_uni(&rnd_uni_init);
	vector<short> neighbor_label;		//存储i节点邻居label	
	double ind_func = namda[0] * ind.y_obj[0] + namda[1] * ind.y_obj[1];

	short neighborsize = node[RandGene].neighbours.size();

	if (neighborsize > 1)
	{
		for (short nei = 0; nei < neighborsize; nei++)
		{
			short neighbornode = node[RandGene].neighbours[nei];
			neighbor_label.push_back(ind.x_var[neighbornode]);
		}

		vector<short>::iterator iter;

		sort(neighbor_label.begin(), neighbor_label.end());
		iter = unique(neighbor_label.begin(), neighbor_label.end());
		neighbor_label.erase(iter, neighbor_label.end());
		//		vector<vector<short> >neighborhood;//neighborhood[0]不能表示一个vector，shit
		//		vector<vector<double> >func;

		double minaggregatefun = 1000000;
		short nodeindex;
		for (i = 0; i < neighbor_label.size(); i++)
		{
			vector<short> neighborhood;
			vector<double> func;
			func.push_back(0.0);
			func.push_back(0.0);
			neighborhood = ind.x_var;
			neighborhood[RandGene] = neighbor_label[i];
			objectives(neighborhood, func);
			double aggregatefun = func[0] * namda[0] + func[1] * namda[1];
			if (minaggregatefun > aggregatefun)
			{
				minaggregatefun = aggregatefun;
				nodeindex = i;
			}
		}
		if (minaggregatefun < ind_func)
		{
			ind.x_var[RandGene] = neighbor_label[nodeindex];
		}

	}//end if neighbor.size > 1
	else {
		if (neighborsize == 1)
		{
			ind.x_var[node[RandGene].neighbours[0]] = ind.x_var[RandGene];
		}
		else {
			cout << "第" << i + 1 << "个节点无连接" << endl;
		}

	}

}


#endif