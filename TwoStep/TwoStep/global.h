#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <windows.h>
#include <iomanip>

using namespace std;

#include "random.h"

#define GNEx_nodes 128
#define GNEx_edge 1024
#define GNEx_QR 1.0/(2*GNEx_edge)

#define INF 1.0e+30
#define archive_size 500   /* set capacity of archive */

extern int     numVariables;
extern int     numObjectives;
extern int     popsize;
extern int     niche;
extern int     mt;
extern int     max_gen;
extern int     Tr;
extern int     nr;

//*
extern char    strFunctionType[256];
extern char    FileName[256];
extern char    FilePath[256];
extern char    NetworkName[256];
extern char    NetworkLabel[256];
extern char    LabelName[256];

// ideal point used in decomposition methods
extern double  *idealpoint;

// parameters for random number generation
extern int     seed;
extern long    rnd_uni_init;

extern int optimization;	   /* set optimization type, 0 for min, 1 for max */
extern double alpha;
extern double r;
extern double delta;
extern double deta;
extern const double perterbation;
extern double mutate_posibility ;
/*--实验发现对于人造数据集必须是小于等于0.01才可以，而实际数据集必须0.1才可以---*/
extern const double cross_posibility;
//extern vector<vector<double> >front; // store removed PF

extern vector<vector<int> >chrom;    // store the corresponding chrom
extern vector<double> NMI;
extern vector<double> pm_NMI;
extern vector<double> pm_modularity;
extern vector<double> Thirty_Run_NMI;
extern vector<double> Thirty_Run_modularity;
extern vector<int> clusters;
extern double NMIavg;
extern double NMImax;
extern double Qmax;
extern double Qavg;
extern int SignedFlag;		//1 for signed networks, 0 for unsigned networks

extern vector<double> modularity;
//extern unsigned int nondomCtr = 0;		     /* number of nondominated solutions in archive      */

struct network 
{
	vector<int> neighbours;		//positive neighbors
	vector<int> neighbours_n;	//negative neighbors
	int degree;					//positive links
	int degree_n;				//negative links
};

//extern int **AdjacentMatrix ;
extern char **AdjacentMatrix   ;
extern network *node ;



#endif
