


/***************************************************************************************

			版本1参考师兄程序，版本2为自己所写的代码(对版本1的优化)

	1. MOEAD: generate an initial population x(上标1,2,3,4,....N);粒子数为N

	2. 被分成N个子问题;indiv 对应着numobjectives?;粒子对应population,种群的indiv对应node,每个population 对应着两个函数值

	3. 每一个x(上标)对应节点个数个x(下标);
	
	3. recombination 头文件中有改进的一段代码;

	3  numvariables为node的个数;pbest.push_back(population[i].indiv),pbest对应indiv,gbest对应population

	4. 为每个种群(例如种群i)基于种群间的权向量之间的距离选择niche个邻居，更新该种群的x_var;x_var 代表位置？

	5. 每个种群对应一个weightvector;

	6. 编程验证a,&a的区别;

	7. 虚函数，this指针的用法;

	8. 容器和数组的区别;

	9. update_reference()函数的位置??????

	10. store_pbest() 使用方法

	11. void TMOEAD::update_pbests_child(TIndividual &ind, int index)

	12.if (strcmp(name[n] , "#")) //这句话不能这么写 if(name[n] != "#")//	
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

xvar随机赋值,邻居加权最小，权重2  100代时 0.3010,200代时 0.31616均匀性也很好，300代 0.298765

运行代数多的话基本上每次运行都会出来很多点，代数少的话有可能哪次运行就出现点数比较少

变异概率越小，得到的点越少，

一个子问题代表一个粒子（n个顶点(xvar)）

population每一行代表一个粒子;每个粒子有no. xvar xvar的维数即是顶点的个数
****************************************************************************************/
/*   MOEA/D DE
void update_problem(TIndividual &indiv, int id)
{
	double random=rnd_uni(&rnd_uni_init);
	double *f1=new double[pops];
	double *f2=new double[pops];
	int c=0;
	if (random<deta)  //p取邻居
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
//更新权重
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
			pop_x.push_back(population[i]); //存储非支配个体
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


交叉变异
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
