


/***************************************************************************************

			�汾1�ο�ʦ�ֳ��򣬰汾2Ϊ�Լ���д�Ĵ���(�԰汾1���Ż�)

	1. MOEAD: generate an initial population x(�ϱ�1,2,3,4,....N);������ΪN

	2. ���ֳ�N��������;indiv ��Ӧ��numobjectives?;���Ӷ�Ӧpopulation,��Ⱥ��indiv��Ӧnode,ÿ��population ��Ӧ����������ֵ

	3. ÿһ��x(�ϱ�)��Ӧ�ڵ������x(�±�);
	
	3. recombination ͷ�ļ����иĽ���һ�δ���;

	3  numvariablesΪnode�ĸ���;pbest.push_back(population[i].indiv),pbest��Ӧindiv,gbest��Ӧpopulation

	4. Ϊÿ����Ⱥ(������Ⱥi)������Ⱥ���Ȩ����֮��ľ���ѡ��niche���ھӣ����¸���Ⱥ��x_var;x_var ����λ�ã�

	5. ÿ����Ⱥ��Ӧһ��weightvector;

	6. �����֤a,&a������;

	7. �麯����thisָ����÷�;

	8. ���������������;

	9. update_reference()������λ��??????

	10. store_pbest() ʹ�÷���

	11. void TMOEAD::update_pbests_child(TIndividual &ind, int index)

	12.if (strcmp(name[n] , "#")) //��仰������ôд if(name[n] != "#")//	
	       name[emptynode[m]] = "#";	//shit ass

double distancevector (vector <double> vec1,vector <double> vec2)
{
	double sum=0;
	int dim=vec1.size();
	for (int i=0;i<dim;i++)
	{
		sum+=(vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
	}
	return sqart(sum);
}

xvar�����ֵ,�ھӼ�Ȩ��С��Ȩ��2  100��ʱ 0.3010,200��ʱ 0.31616������Ҳ�ܺã�300�� 0.298765

���д�����Ļ�������ÿ�����ж�������ܶ�㣬�����ٵĻ��п����Ĵ����оͳ��ֵ����Ƚ���

�������ԽС���õ��ĵ�Խ�٣�

һ�����������һ�����ӣ�n������(xvar)��

populationÿһ�д���һ������;ÿ��������no. xvar xvar��ά�����Ƕ���ĸ���
****************************************************************************************/
/*   MOEA/D DE
void update_problem(TIndividual &indiv, int id)
{
	double random=rnd_uni(&rnd_uni_init);
	double *f1=new double[pops];
	double *f2=new double[pops];
	int c=0;
	if (random<deta)  //pȡ�ھ�
	{
		for (int j=0; j<niche; j++)
		{
			f1[j]=scalar_func(indiv.y_obj, population[j].namda, indivpoint);
			f2[j]=scalar_func(population[j].indiv.y_obj, population[j].namda, indivpoint);
			if (f1[j]<=f2[j])
			c=c+1;
		}
		if (c<=nr)
		{
			for (j=0;j<niche;j++)
			{
				if (f1[j]<f2[j])
				{
					population[j].indiv=indiv;
				}
			}
		}
		else
		{
			double x=new double[niche];
			int idx=new int[niche];
			int idx1=0;
			for (j=0;j<niche;j++)
			{
				if (f1[j]<f2[j])
				{
					x[idx1]=distanceVector(population[j].indiv.y_obj,indiv.y_obj)
					idx[idx1]=j;
					idx1++;
				}
			}
			minfastsort(x,idx,niche,nr);
			for(int k=0; k<nr; k++)   
				{
					//population[id].tabler.push_back(idx[k]);
					population[idx[k]].indiv=indiv;
				}
		}
	}//end if
	else
	{
		for (int j=0; j<pops; j++)
		{
			f1[j]=scalar_func(indiv.y_obj, population[j].namda, indivpoint);
			f2[j]=scalar_func(population[j].indiv.y_obj, population[j].namda, indivpoint);
			if (f1[j]<=f2[j])
			c=c+1;
		}
		if (c<=nr)
		{
			for (j=0;j<pops;j++)
			{
				if (f1[j]<f2[j])
				{
					population[j].indiv=indiv;
				}
			}
		}
		else
		{
			double x=new double[pops];
			int idx=new int[pops];
			int idx1=0;
			for (j=0;j<pops;j++)
			{
				if (f1[j]<f2[j])
				{
					x[idx1]=distanceVector(population[j].indiv.y_obj,indiv.y_obj)
					idx[idx1]=j;
					idx1++;
				}
			}
			minfastsort(x,idx,pops,nr);
			for(int k=0; k<nr; k++)   
				{
					//population[id].tabler.push_back(idx[k]);
					population[idx[k]].indiv=indiv;
				}
		}
	}
}
for(int i=0; i<=sd; i++)
	{
        TSOP sop;
		double a;
		a=1.0*i/sd;
		sop.namda.push_back(a);
		sop.namda.push_back(1-a);
		population.push_back(sop); 		
	}
//����Ȩ��
void update_weightvector(int id)
{
	for (int i=0;i<pops;i++)
	{
		double acc=population[i].namda[0];
		for (int j=1;j<numObjectives;j++)
		{
			acc=acc*population[i].namda[j];
		}
		if (m!=0)
		{
			for (j=0;j<numObjectives;j++)
			{
				double sum=0.0;
				sum=sum+population[i].namda[j];
			}
			for (j=0;j<numObjectives;j++)
			{
				population[i].namda[j]=1.0/population[i].namda[j];
				population[i].namda[j]=population[i].namda[j]/sum;
			}
		}
	}
}

void non_dominated ()
{
	vector <TSOP> pop_x; 
	for (int i=0;i<pops;i++)
	{
		population[i].CNT1=0;
		population[i].n=0;
		for (int j=0;j<pops;j++)
		{
			population[i].sp[j]=0;
		}
	}
	for ( i=0;i<pops;i++)
	{
		for ( j=0; j<pops;j++)
		{
			if (j!=i)
			{
				if (population[i].indiv.y_obj[0]<population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<population[j].indiv.y_obj[1] ||
					population[i].indiv.y_obj[0]<=population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<population[j].indiv.y_obj[1] ||
					population[i].indiv.y_obj[0]<population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<=population[j].indiv.y_obj[1])
				{
					population[i].sp[population[i].CNT1]=j;
					population[i].CNT1++;
				}
				else if (population[i].indiv.y_obj[0]<population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<population[j].indiv.y_obj[1] ||
					population[i].indiv.y_obj[0]<=population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<population[j].indiv.y_obj[1] ||
					population[i].indiv.y_obj[0]<population[j].indiv.y_obj[0]&&population[i].indiv.y_obj[1]<=population[j].indiv.y_obj[1])
				{
					population[i].n++;
				}
		}
	}
	for (i=0;i<pops;i++)
	{
		if (population[i].n==0)
		{
			pop_x.push_back(population[i]); //�洢��֧�����
		}
	}
	fstream fout_front;
	sprintf(saveFilename,"D:\\%d.txt",rutimes+1);
	fout_front.open(saveFilename,ios_base::out);
	for(int n=0; n<pop_x.size(); n++)
	{
		for(int k=0;k<numObjectives;k++)
			fout_front<<pop_x[n][k]<<"  ";
		fout_front<<endl;
	}
	fout_front.close();
}
// sigma method in PSO chose the gbest


�������
template <class T>
void cross_mutation (T &ind1, T &ind2)
{
	for (i=0; i<int(N/2); i++)
	{
		parent1 = int(N*rnd_uni(&rnd_uni_init));
		r0 = rnd_uni(&rnd_uni_init);
		if (r0 <= 0.9)
		{
			parent2 = int(N*rnd_uni(&rnd_uni_init));
			while (parent1 == parent2)
			{
				parent2 = int(N*rnd_uni(&rnd_uni_init));
			}

		}
	}
}

*/
