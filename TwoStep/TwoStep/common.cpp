#include "common.h"

int flag(vector<double> F1,vector<double> F2)
{
	if (F1[0] < F2[0] && F1[1] < F2[1]) return 1;

	if (F1[0] < F2[0] && F1[1] == F2[1]) return 1;

	if (F1[0] == F2[0] && F1[1] < F2[1]) return 1;

	if (F1[0] == F2[0] && F1[1] == F2[1]) return 2;

	return 3;
}



double max_var(vector<double> &x_var)
{
	double temp = x_var[0];
	for (int i=1; i<x_var.size(); i++)
	{
		if (temp < x_var[i])
		{
			temp = x_var[i];
		}
	}
	return temp;
}

double min_var(vector<double> &x_var)
{
	double temp = x_var[0];
	for (int i=1; i<x_var.size(); i++)
	{
		if (temp > x_var[i])
		{
			temp = x_var[i];
		}
	}
	return temp;
}

double compute_tr(int iter, int n, double r)
{
	double tr;
	double tr1;

	tr1=(-20)*(iter/max_gen-r);
	tr=0.4*n/(1+exp(tr1));
	tr=int (tr);

	return tr;
}
double distanceArray(double vec1[], double vec2[], int dim)
{
    double sum = 0;

	for(int n=0; n<dim; n++)
	{
		sum+= (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	}

	return sqrt(sum);
}

double distanceVector(vector <double> &vec1, vector <double> &vec2)
{
    double sum = 0;

	for(int n=0; n<vec1.size(); n++)
	{
		sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	}

	return sqrt(sum);
}


double norm_vector(vector <double> &x)
{
	double sum = 0;

	for(int i=0;i<x.size();i++)
	{
		sum = sum + x[i]*x[i];
	}

    return sqrt(sum);
}

double sum_vector(vector<double>&vec)
{
	double sum = 0;

	for(int i=0;i<vec.size();i++)
	{
		sum = sum + vec[i];
	}

    return sum;
}

double innerproduct(vector <double>&vec1, vector <double>&vec2)
{
    double sum = 0;

	for(int i=0; i<vec1.size(); i++)
	{
		sum+= vec1[i]*vec2[i];
	}

	return sum;
}

void minfastsort(double x[], int idx[], int n, int m)
//void minfastsort(vector <double> x, vector <int> idx, int n, int m)
{
    for(int i=0; i<m; i++)
	{
		for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
				double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}

void ReadFile(char *filepath,char** p,int row,int col)

{
	ifstream input(filepath,ios::in);
	
	if(input.fail())
	{
		cout<<"Oops!Couldn't open data file\n";
		exit(1);
	}
	
	for (int i = 0;i < row;i++ )
	{
		for (int j = 0;j < col;j++ )
		{
			input >> p[i][j];
			input.get();
		}
	}
	
	input.close();
}

void labpos(vector<int> gene,network *node)
{
	for (int i=0;i<numVariables;i++)
	{
		int NeighborSize = node[i].neighbours.size();

		if (NeighborSize == 0)
		{
			gene[i] = i;
		} 
		else
		{
			int temp;
			temp = int (rnd_uni(&rnd_uni_init)*NeighborSize);
			gene[i] = gene[node[i].neighbours[temp]];
		}
	}
}

void NodeInformation()
{
	for ( int i = 0; i<numVariables; i++ )
	{
		int positive_sum = 0;
		int negative_sum = 0;
		
		for ( int j = 0;j < numVariables;j++ )
		{
// 			if ( AdjacentMatrix[i][j] == 1 )
// 			{
// 				node[i].neighbours.push_back(j);
// 				positive_sum += AdjacentMatrix[i][j];
// 			}
// 			else if (AdjacentMatrix[i][j] == -1)
// 			{
// 				//abs for int,fabs for float and labs for long
// 				node[i].neighbours_n.push_back(j);
// 				negative_sum += abs(AdjacentMatrix[i][j]);
// 			}	
			

			if (AdjacentMatrix[i][j]==0x31)
			{
				node[i].neighbours.push_back(j);
				positive_sum += 1;
			}
			else 
			{
				if (AdjacentMatrix[i][j]==0x2d)	
				{
					node[i].neighbours_n.push_back(j);
					negative_sum += 1;
				}
			}

		}
		node[i].degree = positive_sum;
		node[i].degree_n = negative_sum;
	}
}


double calc_NMI(vector<int> gene,char *filepath)

{
	int i,j,k,m;
	double Temp_NMI = 0.0;
	int* p = new int[numVariables];

	ifstream input(filepath,ios::in);   //open the content in the filepath 
		
	if(input.fail())
	{
		cout<<"Oops!Couldn't open real label file\n";
	//	exit(1);
	}
		
	for ( i = 0;i < numVariables;i++ )
	{
		input >> p[i];           //input the content of the input in the p matrix
		input.get();
	}

	vector<int> copylabel;
	copylabel = gene;

	vector<vector<int> >A;  //the detected community
	vector<vector<int> >B;  //the real community

	for ( i = 0;i < numVariables;i++ )
	{	//calculate A				
		if ( copylabel[i] != -1 )
		{
			vector<int> s1;	//store every component in each cluster
			s1.push_back(i);
			for ( j = i+1;j < numVariables;j++ )
			{
				if ( copylabel[i] == copylabel[j] )
				{
					s1.push_back(j);
					copylabel[j] = -1;
				}
			}//end j
			copylabel[i] = -1;
			A.push_back(s1);
		}//end if
			//calculate B
		if ( p[i] != -1 )
		{
			vector<int> s2;	//store every component in each cluster
			s2.push_back(i);
			for ( k = i+1;k < numVariables;k++ )
			{
				if ( p[i] == p[k] )
				{
					s2.push_back(k);
					p[k] = -1;
				}
			}//end k
			p[i] = -1;
			B.push_back(s2);
		}//end if
	}//end i		
	//cout<<"分成了"<<A.size()<<"类";
/*/////////////////////////////////////////
		cout<<"detected A"<<endl;
	for ( i = 0;i < A.size();i++ )
	{
		cout<<"第 "<<i+1<<" 类:"<<endl;
		for ( k = 0;k < A[i].size();k++ )
		{
		//	cout<<A[i][k]+1<<' ';
			cout<<A[i][k]<<' ';
		}
		cout<<endl;
	}cout<<endl;
/*//////////////////////////////////////////*/

	/*/////////////////////////////////////////
		cout<<"real ones"<<endl;
	for ( i = 0;i < B.size();i++ )
	{
		cout<<"第 "<<i+1<<" 类:"<<endl;
		for ( k = 0;k < B[i].size();k++ )
		{
		//	cout<<A[i][k]+1<<' ';
			cout<<B[i][k]<<' ';
		}
		cout<<endl;
	}cout<<endl;
/*//////////////////////////////////////////*/
	double x = 0.0;
	double y = 0.0;
	int *ci = new int[A.size()];//sum of row
	int *cj = new int[B.size()];//sum of col
//	vector<vector<int> >C;		//the fuzzy matrix
	vector<vector<int> >C(A.size(),vector<int>(B.size()));

	for ( i = 0;i < A.size();i++ )
		ci[i] = 0;
	for ( j = 0;j < B.size();j++ )
		cj[j] = 0;

	for ( i = 0;i < A.size();i++ )
		for ( j = 0;j < B.size();j++ )
			C[i][j] = 0;

	for ( i = 0;i < A.size();i++ )
	{
		for ( j = 0;j < B.size();j++ )
		{
			for ( k = 0;k < A[i].size();k++ )
			{
				for ( m = 0;m < B[j].size();m++ )
				{
					if ( A[i][k] == B[j][m] )
					{
						C[i][j] += 1;
						break;
					}//end if
				}//end m
			}//end k
		}//end j	
	}//end i

	for ( i = 0;i < A.size();i++ )
	{
		for ( j = 0;j < B.size();j++ )
		{
			ci[i] += C[i][j];
		}
	}

	for ( j = 0;j < B.size();j++ )
	{
		for ( k = 0;k < A.size();k++ )
		{
			cj[j] += C[k][j];
		}
	}
	
	double temp1 = 0.0;

	for ( i = 0;i < A.size();i++ )
	{
		temp1 += 1.0*ci[i]*log10(1.0*ci[i]/numVariables);
		for ( j = 0;j < B.size();j++ )
		{
			if ( C[i][j] != 0 )
			{
				x += 1.0*C[i][j]*log10(1.0*C[i][j]*numVariables/(ci[i]*cj[j]));
			}
		}
		
	}

	for ( j = 0;j < B.size();j++ )
	{
		y += 1.0*cj[j]*log10(1.0*cj[j]/numVariables);
	}

	Temp_NMI = -2*x/(temp1 + y);
	

	delete[] p;
	delete[] ci;
	delete[] cj;
	A.clear();
	B.clear();
	C.clear();

	return Temp_NMI;
}

bool check_label()
{
	if (strlen(NetworkLabel) > 0)
		return true;
	else return false;
}

//template<class T>
void rehearse(vector<int> x)
{
	vector<int> copylabel;

	copylabel = x;
	
	vector<vector<int> >A;  //the detected community
	
	for ( int i = 0;i < numVariables; i++ )
	{	//calculate A				
		if ( copylabel[i] != -1 )
		{
			vector<int> s1;	//store every component in each cluster
			s1.push_back(i);
			for ( int j = i+1;j < numVariables;j++ )
			{
				if ( copylabel[i] == copylabel[j] )
				{
					s1.push_back(j);
					copylabel[j] = -1;
				}
			}//end j
			copylabel[i] = -1;
			A.push_back(s1);
		}//end if
	}//end i	
	
	for (int t_row = 0; t_row < A.size(); t_row++)
	{
		for (int t_col = 0; t_col < A[t_row].size(); t_col++)
		{
			x[A[t_row][t_col]] = t_row+1;
		}
	}
	A.clear();
	copylabel.clear();
}

/*----print the detected partitions for----
  ----the unknown structure networks------*/
void PrintPartition(int* gene,char *filepath)
{
	int i,j,k;

	int* copylabel = new int[numVariables];

	for( j = 0;j < numVariables;j++ )
		copylabel[j] = gene[j];

// 	for (int emptynode = 0; emptynode < nodes; emptynode++)
// 	{
// 		if (node[emptynode].degree == 0)
// 			copylabel[emptynode] = 0;
// 	}

	vector<vector<int> >A;  //the detected community

	for ( i = 0;i < numVariables;i++ ){	//calculate A				
		if ( copylabel[i] != -1 ){
			vector<int> s1;	//store every component in each cluster
			s1.push_back(i);
			for ( j = i+1;j < numVariables;j++ ){
				if ( copylabel[i] == copylabel[j] ){
					s1.push_back(j);
					copylabel[j] = -1;
				}
			}//end j
			copylabel[i] = -1;
			A.push_back(s1);
		}//end if
	}//end i	
	
	delete [] copylabel;
	
	int* p = new int[numVariables];
	
	for ( i = 0;i < A.size();i++ ){
		for ( k = 0;k < A[i].size();k++ ){
			p[A[i][k]] = i + 1;
		}
	}
	
//	cout<<"检测出 "<<A.size()<<" 类社区"<<endl;
	
	ofstream output(filepath,ios::out);
	output<<"*vertices "<<numVariables<<endl;
				
	for ( i = 0;i < numVariables;i++ )
		output<<p[i]<<endl;
	
	delete [] p;
	output.close();
	A.clear();
}



int deltafun(int a, int b)
{
	if (a == b)
		return 1;
	else return 0;
}


double calcQ(vector<int> pos, int &Anum1, int &countn, int &countp)
{
	double R;
	unsigned int i, j, k, m;
	//int SignedFlag = 1;
	int* copylabel = new int[numVariables];
	
	for( j = 0;j < numVariables;j++ )
		copylabel[j] = pos[j];
	
	vector<vector<int> >A;  //the detected community
	
	for ( i = 0;i < numVariables;i++ )
	{	//calculate A				
		if ( copylabel[i] != -1 )
		{
			vector<int> s1;	//store every component in each cluster
			s1.push_back(i);
			for ( j = i+1;j < numVariables;j++ )
			{
				if ( copylabel[i] == copylabel[j] )
				{
					s1.push_back(j);
					copylabel[j] = -1;
				}
			}//end j
			copylabel[i] = -1;
			A.push_back(s1);
		}//end if
	}
	Anum1 = A.size();
//	cout<<"检测出 "<<A.size()<< "类";

	for (i=0; i<A.size(); i++)
	{
	//	cout << "第" << i <<"类"<< " ";
		for (j=0; j<A[i].size(); j++)
		{
		//	cout << "第" <<i<<"类中元素个数=" << A[i].size() << endl;
			for (k=i+1; k<A.size(); k++)
			{
				for (m=0; m<A[k].size(); m++)
				{
					int a1 = A[i][j];
					int b1 = A[k][m];
			//		cout << "第" << i << "类中的第" << j << "个元素和第" << k << "类中的第" << m << "个元素比较" <<endl;
// 					if (AdjacentMatrix[a1][b1] == 1)
// 					{
// 						countn++;
// 					}

					if (AdjacentMatrix[a1][b1]==0x31)
					{
						countn++;
					}

				}
			}
		}
	}
	
	for (i=0; i<A.size(); i++)
	{
		for (j=0; j<A[i].size(); j++)
		{
			for (k=j+1; k<A[i].size(); k++)
			{
				int a2 = A[i][j];
				int b2 = A[i][k];
		//		cout << "第" << i << "类中的第" << j << "个元素和第" << i << "类中的第" << k << "个元素比较" <<endl;
// 				if (AdjacentMatrix[a2][b2] == -1)
// 				{
// 					countp++;
// 				}

				if (AdjacentMatrix[a2][b2]==0x2d)
				{
					countp++;
				}

			}
		}
	}
	
/*	fstream fout_front;
//	sprintf(saveFilename,"F:\\电脑资料\\hs\\SP\\%d.txt",num+1);
	fout_front.open(saveFilename,ios_base::out);
	for(int n=0; n<numVariables; n++)
	{
			fout_front<<pos[n]<<"  ";
	}
	fout_front.close();*/


	A.clear();
	delete[] copylabel;

	/*if (SignedFlag == 0)
	{
		if(!strcmp(strReal,"FB50.txt")) R = FB50_QR;
		if(!strcmp(strReal,"karate.txt")) R = karate_QR;
		if(!strcmp(strReal,"dolphin.txt")) R = dolphin_QR;
		if(!strcmp(strReal,"SFI.txt")) R = SFI_QR;
		if(!strcmp(strReal,"football.txt")) R = football_QR;
		if(!strcmp(strReal,"netscience.txt")) R = netscience_QR;
		if(!strcmp(strReal,"netscience_remove.txt")) R = netscience_remove_QR;
		if(!strcmp(strReal,"power.txt")) R = power_QR;
		if(!strcmp(strReal,"power_remove.txt")) R = power_remove_QR;
		if(!strcmp(strReal,"hepth_adj.txt")) R = hepth_QR;	
		if(!strcmp(strReal,"PGP_adj.txt")) R = PGP_QR;
		if(!strcmp(strReal,"G.txt")) R = MIDV30_QR;
		if(!strcmp(strReal,"MID3.0_remove.txt")) R = MIDV30_remove_QR;

		
		double Q = 0;
		for ( i = 0;i < numVariables;i++ ){
			int di = node[i].degree;
			int u = pos[i];
			for ( j = 0;j < numVariables;j++ ){
				int dj = node[j].degree;
				int v = pos[j];
				Q += (AdjacentMatrix[i][j] - di * dj * R)*deltafun(u,v);
			}
		}
		Q *=  R;
		return Q;
	} */
	//else if (SignedFlag == 1)
	if (SignedFlag == 1)
	{		
		double Q = 0;
		double positive_degree = 0, negative_degree = 0;
		for (i = 0;i < numVariables;i++)
		{
			positive_degree += node[i].degree;
			negative_degree += node[i].degree_n;
		}
		R = positive_degree + negative_degree;
		R = 1.0/R;
		positive_degree = 1.0/positive_degree;
		negative_degree = 1.0/negative_degree;
		for ( i = 0;i < numVariables;i++ ){
			int idp = node[i].degree;
			int idn = node[i].degree_n;
			int u = pos[i];
			for ( j = 0;j < numVariables;j++ ){
				int jdp = node[j].degree;
				int jdn = node[j].degree_n;
				int v = pos[j];
// 				Q += (AdjacentMatrix[i][j]-(idp*jdp*positive_degree-idn*jdn*negative_degree))*deltafun(u,v);

				int aij = 0;
				if (AdjacentMatrix[i][j]==0x2d)
				{
					aij = -1;
				}
				else
				{
					if (AdjacentMatrix[i][j]==0x31)
					{
						aij = 1;
					}
				}

				Q += (aij-(idp*jdp*positive_degree-idn*jdn*negative_degree))*deltafun(u,v);

			}
		}
		Q *=  R;
	return Q;
	}
	return -1;
}
/*--only for two objectives scinorio--*/
/*-----remove the dominated solutions----*/
int donimate_judge(vector<double> pf1, vector<double> pf2)
{
	//return 1: point pf2 donimates pf1;return -1: point pf1 donimates pf2
	int biger = 0;
	int smaller = 0;
	for (int i=0; i<numObjectives; i++)
	{
		if (pf1[i]>=pf2[i])	biger++;
		if (pf1[i]<=pf2[i]) smaller++;		
	}
	if (smaller == numObjectives)	return -1;
	else if (biger == numObjectives) return 1;
	else return 0;
}

void GetNondominatedSolution(vector<vector <double> > population, vector<vector <double> > &NondomiSolutions, vector<int> &index) 
{
	// population: 表示还没有去掉被支配解的一组解
	int temp,i,j;
	for (i = 0; i < NondomiSolutions.size(); i++)		
		NondomiSolutions[i].clear();

	NondomiSolutions.clear();

	int *num_donimated = new int[popsize];
	for (i = 0; i < popsize; i++)
		num_donimated[i] = 0;

	for (i = 0; i < population.size(); i++)
	{
		 for (j = i + 1; j < population.size(); j++)
		 {
			 temp = donimate_judge(population[i], population[j]);
			 if (temp == -1) //i donimate j
				 num_donimated[j]++;
			 else if (temp == 1)//j donimate i
				 num_donimated[i]++;
		 }

		 if (num_donimated[i] == 0) 
		 {
			 NondomiSolutions.push_back(population[i]);
			 index.push_back(i);
		 }
	}
}
/*--only for two objectives scinorio--*/
/*-----remove same point in the PF----*/
vector<vector<double> > remove_same_point(vector<vector<double> >pop)
{
	vector<double> pop_col1;
	vector<double> pop_col2;
	int i,j;
	for ( i = 0; i < pop.size(); i++)
	{
		pop_col1.push_back(pop[i][0]);
		pop_col2.push_back(pop[i][1]);
	}
	
	int flag[100];
	for (i=0; i<100; i++)
	{
		flag[i] = 1;
	}

	for (i=0; i<pop.size(); i++)
	{
		if (flag[i] != -1)
		{
			for (j=i+1; j<pop.size(); j++)
			{
				if ((pop_col1[i] == pop_col1[j])&&(pop_col2[i] == pop_col2[j]))
				{
					flag[j] = -1;
				}
			}
		}
	}

	vector<vector<double> >frontRemove;
	vector<int> index;
	for (i=0; i<100; i++)
	{
		if (flag[i] == 1)
		{
			vector<double> temp;
			temp.push_back(pop_col1[i]);
			temp.push_back(pop_col2[i]);
			frontRemove.push_back(temp);
			index.push_back(i);
		}
	}
	return frontRemove;
}
/*	vector<double>::iterator pop_col1_iter;
	vector<double>::iterator pop_col2_iter;

	sort(pop_col1.begin(),pop_col1.end());
	pop_col1_iter = unique(pop_col1.begin(),pop_col1.end());
	pop_col1.erase(pop_col1_iter,pop_col1.end());

	sort(pop_col2.begin(),pop_col2.end());
	pop_col2_iter = unique(pop_col2.begin(),pop_col2.end());
	pop_col2.erase(pop_col2_iter,pop_col2.end());

	if (pop_col2.size() != pop_col1.size())
	{
		cerr<<"PF面存在支配解"<<endl;
	//	exit(1);
	}

	reverse(pop_col2.begin(),pop_col2.end());
	
	vector<vector<double> >front;
	for (int j = 0; j < pop_col1.size(); j++)
	{
		vector<double> temp;
		temp.push_back(pop_col1[j]);
		temp.push_back(pop_col2[j]);

		front.push_back(temp);
	}

	return front;
}*/

void perturbation(int a[],int N)
{
    int i;
    int index;
    srand((unsigned)time(NULL));
    for (i = N - 1; i > 0; i--)
    {
        index = rand() % i;
		swap(a[i],a[index]);
    }
}

void random_init(vector <int> gene)
{
	int i;
	for (i=0; i<numVariables; i++)
	{
		gene[i]=int (rnd_uni(&rnd_uni_init)*i); 
	}
}

void IGLP(vector<int> gene, network* node)
{   
	int n,i,j,k,l;

	for (i=0; i<numVariables; i++)
	{
		gene[i]=i;
	}

	for (n=0; n<5; n++)
	//for ( n = 0;n < 1;n++ )
	{
		for (i=0; i<numVariables; i++)
		{	
			int NeighborSize = node[i].neighbours.size();
			if (NeighborSize == 0)
			{
//				gene[i] = 0;
			} 
			else
			{
				if (NeighborSize == 1)
				{
					gene[i] = gene[node[i].neighbours[0]];
				//	gene[node[i].neighbours[0]] = gene[i]; //这个实验结果没有上面那个好
				} 
				else
				{
					int	sum = 0;
					int maxr = -1;//record index of i's neighbour which ...	
					int label = -1;
					int	temp = 1;
					
					for (j=0; j<NeighborSize; j++)
					{
						int counter = 1;//record no. of nodes that has same label with j
						
						for (k=j+1; k<NeighborSize; k++)
						{
							int p = gene[node[i].neighbours[j]];
							int q = gene[node[i].neighbours[k]];
							if(p == q)  
							{
								counter++;
							}
						}//end k
						
						if (temp < counter)
						{
							maxr = j;
							temp = counter;
						}
					}//end j
					
					for (l=0; l<NeighborSize; l++)
					{
						int u = gene[node[i].neighbours[l]];
						int v = gene[i];
						if (u == v)
						{
							label = u;
						}
					}//end l
					if ((label != -1) && (maxr == -1))
					{
						gene[i] = label;
					}
					else 
					{
						if (maxr != -1)
						{
							gene[i] = gene[node[i].neighbours[maxr]];
						} 
						else
						{
							int randneighbor = (int)(rnd_uni(&rnd_uni_init)*NeighborSize);
							gene[i] = gene[node[i].neighbours[randneighbor]];
						}
					}
				}
			}
		}//end i
	}//end n
}


/*
 *	如果邻居有大多数（最少两个）个体有相同标号i，则
 *  将该节点的类标换成i，如果没有这样的标号则看有没有
 *  邻居节点标号和该节点标号相同，有就选择它，否则就
 *  随机选取一个邻居的编号
 */
void IGLP1(vector<int> gene, network* node)

{   
	int n, i, j, k, l;
	int *a = new int[numVariables];
	for (int a_row = 0; a_row < numVariables; a_row++)
		a[a_row] = a_row;
	
	perturbation(a, numVariables);

	for ( n = 0;n < 5;n++ )
	{
		for ( i = 0;i < numVariables;i++ )
		{	
			int randnode = a[i];
			int	sum = 0;
			int maxr = -1;//record index of i's neighbour which ...	
			int label = -1;
			int	temp = 1;
			
//			if(node[i].neighbours.size() > 1)
			if(node[randnode].neighbours.size() > 1)
			{
//				for ( j = 0;j < node[i].neighbours.size() - 1;j++ )
				for ( j = 0;j < node[randnode].neighbours.size() - 1;j++ )
				{
					int counter = 1;//record no. of nodes that has same label with j
					
//					for ( k = j + 1;k < node[i].neighbours.size();k++ )
					for ( k = j + 1;k < node[randnode].neighbours.size();k++ )
					{
// 						int p = gene[node[i].neighbours[j]];
// 						int q = gene[node[i].neighbours[k]];
						int p = gene[node[randnode].neighbours[j]];
 						int q = gene[node[randnode].neighbours[k]];
						if( p == q )  counter++;
					}//end k
					
					if ( temp < counter )
					{
						maxr = j;
						temp = counter;
					}
				}//end j
				
//				for ( l = 0;l < node[i].neighbours.size();l++ )
				for ( l = 0;l < node[randnode].neighbours.size();l++ )
				{
//					int u = gene[node[i].neighbours[l]];
					int u = gene[node[randnode].neighbours[l]];
//					int v = gene[i];
					int v = gene[randnode];
					if ( u == v )
					{
						label = u;
					}
				}//end l
				if (label != -1 && maxr == -1)
				{
//					gene[i] = label;
					gene[randnode] = label;
				}
				else 
				{
					if (maxr != -1)
					{
//						gene[i] = gene[node[i].neighbours[maxr]];
						gene[randnode] = gene[node[randnode].neighbours[maxr]];
					} 
					else
					{
//						int randneighbor = rnd(0,1)*node[i].neighbours.size();
//						gene[i] = gene[node[i].neighbours[randneighbor]];
						int randneighbor = rnd_uni(&rnd_uni_init)*node[randnode].neighbours.size();
						gene[randnode] = gene[node[randnode].neighbours[randneighbor]];
					}
				}				
			}//end if
			else 
			{
// 				if (node[i].neighbours.size() == 1)
// 				{
// 					gene[i] = gene[node[i].neighbours[0]];
// 				} 
				if (node[randnode].neighbours.size() == 1)
				{
						gene[randnode] = gene[node[randnode].neighbours[0]];
 				}
				else
				{
//					cout<<"第"<<i+1<<"个节点无连接"<<endl;
				}
				
			}
		}//end i
	}//end n
	delete [] a;
}


void heuristic_initial(vector<int> gene,network* node)
{
	int thita = 0.2*numVariables;
	for ( int m = 0;m < thita;m++ )
	{
		int position  = int(numVariables * rnd_uni(&rnd_uni_init));
		if (position == numVariables) position -= 1;
		
		if (node[position].neighbours.size() != 0)
		{
			for ( int k = 0;k < node[position].neighbours.size();k++ )
			{			
				gene[node[position].neighbours[k]] = gene[position];
			}
		} 
	}
}

void IGLP2(vector<int> gene,network* node)
{
	int n,i,j,k,l;
	for (i=0;i<numVariables;i++)
	{
		for (j=0;j<numVariables;j++)
		{
		//	if (AdjacentMatrix[i][j]==1)
			if (AdjacentMatrix[i][j]==0x31)
				node[i].neighbours.push_back(j);
		}
	}
	
	for ( n = 0;n < 5;n++ )
	{
		for ( i = 0;i < numVariables;i++ )
		{	
			int NeighborSize = node[i].neighbours.size();
			if (NeighborSize == 0)
			{
//				gene[i] = 0;
				gene[i] = i;
			} 
			else
			{
				if (NeighborSize == 1)
				{
					gene[i] = gene[node[i].neighbours[0]];
				//	gene[node[i].neighbours[0]] = gene[i]; //这个实验结果没有上面那个好
				} 
				else
				{
					int	sum = 0;
					int maxr = -1;//record index of i's neighbour which ...	
					int label = -1;
					int	temp = 1;
					
					for ( j = 0;j < NeighborSize;j++ )
					{
						int counter = 1;//record no. of nodes that has same label with j
						
						for ( k = j + 1;k < NeighborSize;k++ )
						{
							int p = gene[node[i].neighbours[j]];
							int q = gene[node[i].neighbours[k]];
							if( p == q )  counter++;
						}//end k
						
						if ( temp < counter )
						{
							maxr = j;
							temp = counter;
						}
					}//end j
					
					for ( l = 0;l < NeighborSize;l++ )
					{
						int u = gene[node[i].neighbours[l]];
						int v = gene[i];
						if ( u == v )
						{
							label = u;
						}
					}//end l
					if (label != -1 && maxr == -1)
					{
						gene[i] = label;
					}
					else 
					{
						if (maxr != -1)
						{
							gene[i] = gene[node[i].neighbours[maxr]];
						} 
						else
						{
							int randneighbor = (int)(rnd_uni(&rnd_uni_init)*NeighborSize);
							gene[i] = gene[node[i].neighbours[randneighbor]];
						}
					}
				}
			}
		}//end i
	}//end n
}

void IGLP21(vector<int> gene,network* node)
{
	int n,i,j,k,l;
	for (i=0;i<numVariables;i++)
	{
		for (j=0;j<numVariables;j++)
		{
		//	if (AdjacentMatrix[i][j]==1)
			if (AdjacentMatrix[i][j]==0x31)
			{
				node[i].neighbours.push_back(j);
			}
		}
	}
	for (i=0;i<numVariables;i++)
	{
		for (n=0;n<5;n++)
		{
			for (k=0;k<numVariables;k++)
			{
				if (node[i].neighbours.size()>1)
				{
					for (j=0;j<node[i].neighbours.size();j++)
					{
						int counter = 1;//record no. of nodes that has same label with j
						int temp=1;
						int maxr=-1;
						for ( l = j + 1;l < node[i].neighbours.size();l++ )
						{
							int p = gene[node[i].neighbours[j]];
							int q = gene[node[i].neighbours[l]];
							if( p == q )  counter++;
						}//end k
						
						if ( temp < counter )
						{
							maxr = j;
							temp = counter;
						}
					}
					gene[i]=gene[node[j].neighbours[j]];
				}
				else
				{
					gene[i]=gene[node[i].neighbours[0]];
				}
			}
		}
	}
}


