#ifndef __COMMON_H_
#define __COMMON_H_

#include <algorithm>
#include "global.h"

double max_var(vector<int> &x_var);

double min_var(vector<int> &x_var);

double compute_tr(int iter, int n, double r);

double distanceArray(double vec1[], double vec2[], int dim);

double distanceVector(vector <double> &vec1, vector <double> &vec2);

double norm_vector(vector <double> &x);

double sum_vector(vector<double>&vec);

double innerproduct(vector <double>&vec1, vector <double>&vec2);

void minfastsort(double x[], int idx[], int n, int m);

//void minfastsort(vector <double> &x, vector <int> &idx, int n, int m);

//void ReadFile(char *filepath,int** p,int row,int col);

void ReadFile(char *filepath,char** p,int row,int col);

void NodeInformation();

void labpos(vector<int> gene,network *node);

double calc_NMI(vector<int> gene,char *filepath);

bool check_label();

/*----print the detected partitions for----
  ----the unknown structure networks------*/
void PrintPartition(int* gene,char *filepath);

int deltafun(int a, int b);

double calcQ(vector<int> pos, int &Anum1, int &countn, int &countp);

int donimate_judge(vector<double> pf1, vector<double> pf2);

void GetNondominatedSolution(vector<vector <double> > population, vector<vector <double> > &NondomiSolutions, vector<int> &index);

/*--only for two objectives scinorio--*/
/*-----remove same point in the PF----*/
vector< vector<double> > remove_same_point(vector< vector<double> >pop);

void perturbation(int a[],int N);

void random_init(vector <int> gene);

void IGLP(vector<int> gene,network* node);

/*
 *	如果邻居有大多数（最少两个）个体有相同标号i，则
 *  将该节点的类标换成i，如果没有这样的标号则看有没有
 *  邻居节点标号和该节点标号相同，有就选择它，否则就
 *  随机选取一个邻居的编号
 */
void IGLP1(vector<int> gene,network* node);

void heuristic_initial(vector<int> gene,network* node);

void IGLP2(vector<int> gene,network* node);
 
void IGLP21(vector<int> gene,network* node);

#endif

